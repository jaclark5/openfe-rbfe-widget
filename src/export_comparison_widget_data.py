"""Build widget JSON and HTML for OpenFE benchmark comparisons.

Main function: ``export_comparison_widget``.

Arguments
---------
config : str
    YAML config path. Required.

Keyword Arguments
-----------------
output_dir : str
    Output directory.
num_bootstraps : int
    Bootstrap count for summary statistics.
ecdf_bootstraps : int
    Bootstrap count for ECDF confidence intervals.
no_ecdf_ci : bool
    Disable ECDF confidence intervals.
no_structures : bool
    Skip 2D ligand rendering.
data_only : bool
    Write JSON only.
seed : int
    Bootstrap random seed.
dump_counts : bool
    Print verification counts.
reference_model : str | None
    Label used for significance-based bar hatching.

Outputs
-------
comparison_widget_data.json
    Exported widget data.
comparison_widget.html
    Self-contained widget HTML. Omitted when ``data_only=True``.
"""

import base64
import json
import logging
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import urllib.request
import urllib.error

import xml.etree.ElementTree as ET
    
import click
import numpy as np
import pandas as pd
import yaml
from rdkit.Chem import Draw

from cinnabar.compare import compare_and_rank_results
from cinnabar import FEMap, ReferenceState
import openfe
from openfe.utils.atommapping_network_plotting import draw_one_molecule_mapping

from openfe_benchmarks.data import get_data_by_system_name
from openfe_benchmarks.results._benchmark_results import get_benchmark_results

# Setup logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)


def _compute_ecdf_with_bootstrap(
    errors: np.ndarray,
    num_bootstraps: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute ECDF with bootstrap confidence intervals.

    Follows cinnabar.plotting.ecdf_plot's bootstrap algorithm:
    resample errors with replacement, compute ECDF via searchsorted,
    then take 2.5/97.5 percentiles across bootstraps.

    Parameters
    ----------
    errors : np.ndarray
        Absolute errors |calc - exp| for one key, shape (n_edges,).
    num_bootstraps : int
        Number of bootstrap samples for CI bands. 0 disables CI.

    Returns
    -------
    x : np.ndarray
        Sorted error values (step curve x-coordinates).
    y : np.ndarray
        ECDF values i/n (step curve y-coordinates).
    ci_lower : np.ndarray or None
        Lower 2.5th percentile band, or None if num_bootstraps == 0.
    ci_upper : np.ndarray or None
        Upper 97.5th percentile band, or None if num_bootstraps == 0.
    """
    # ECDF: sorted errors and cumulative fractions
    x = np.sort(errors)
    n = len(x)
    y = np.arange(1, n + 1) / n

    if num_bootstraps == 0:
        return x, y, None, None

    # Bootstrap CI following cinnabar.plotting.ecdf_plot
    # Resample errors, compute ECDF at original x via searchsorted
    bootstrap_ecdfs = []
    for _ in range(num_bootstraps):
        sample = np.random.choice(errors, size=n, replace=True)
        sample_sorted = np.sort(sample)
        # ECDF at each x[i]: fraction of sample <= x[i]
        ecdf_at_x = np.searchsorted(sample_sorted, x, side="right") / n
        bootstrap_ecdfs.append(ecdf_at_x)

    bootstrap_ecdfs = np.array(bootstrap_ecdfs)  # shape (num_bootstraps, n)

    # 2.5th and 97.5th percentiles across bootstraps
    ci_lower = np.percentile(bootstrap_ecdfs, 2.5, axis=0)
    ci_upper = np.percentile(bootstrap_ecdfs, 97.5, axis=0)

    return x, y, ci_lower, ci_upper


def _compute_dg_rmse_with_bootstrap(
    dg_centralized_by_key: dict[str, np.ndarray],
    dg_exp: np.ndarray,
    num_bootstraps: int,
) -> dict[str, dict[str, float]]:
    """Compute ΔG RMSE from centralized values with bootstrap CI.

    This is the one hand-rolled statistic (D3 exception). It uses
    cinnabar.compare's joint pairwise resampling (resample
    (prediction, experiment) pairs jointly across models), not
    cinnabar.plotting.ecdf_plot's flat-error-array resampling.

    Parameters
    ----------
    dg_centralized_by_key : dict[str, np.ndarray]
        Centralized computational ΔG per key, shape (n_ligands,).
    dg_exp : np.ndarray
        Experimental ΔG, shape (n_ligands,).
    num_bootstraps : int
        Number of bootstrap samples.

    Returns
    -------
    dict[str, dict[str, float]]
        {key: {"value": rmse, "ci_lower": float, "ci_upper": float}}
    """
    keys = list(dg_centralized_by_key.keys())
    n = len(dg_exp)

    # Compute point estimates
    rmse_values = {}
    for key in keys:
        dg_calc = dg_centralized_by_key[key]
        rmse = np.sqrt(np.mean((dg_calc - dg_exp) ** 2))
        rmse_values[key] = rmse

    # Bootstrap: resample (prediction, experiment) pairs jointly
    bootstrap_rmses = {key: [] for key in keys}

    for _ in range(num_bootstraps):
        indices = np.random.choice(n, size=n, replace=True)
        dg_exp_sample = dg_exp[indices]

        for key in keys:
            dg_calc_sample = dg_centralized_by_key[key][indices]
            rmse_sample = np.sqrt(np.mean((dg_calc_sample - dg_exp_sample) ** 2))
            bootstrap_rmses[key].append(rmse_sample)

    # Compute 2.5th and 97.5th percentiles
    result = {}
    for key in keys:
        ci_lower = np.percentile(bootstrap_rmses[key], 2.5)
        ci_upper = np.percentile(bootstrap_rmses[key], 97.5)
        result[key] = {
            "value": float(rmse_values[key]),
            "ci_lower": float(ci_lower),
            "ci_upper": float(ci_upper),
        }

    return result


def _apply_dg_centralizing_shift(
    dg_computational: pd.Series,
    dg_experimental: pd.Series,
) -> pd.Series:
    """Apply cinnabar.plotting.plot_DGs centralizing shift to ΔG values.

    This applies the same centralizing transformation as cinnabar's plot_DGs
    with centralizing=True: x - mean(x) + experimental_mean, so computational
    and experimental values are on the same offset.

    Parameters
    ----------
    dg_computational : pd.Series
        Computational DG values (kcal/mol), indexed by ligand.
    dg_experimental : pd.Series
        Experimental DG values (kcal/mol), indexed by ligand (same index as computational).

    Returns
    -------
    pd.Series
        Centralized computational DG values.
    """
    # Centralizing shift: x - mean(x) + mean(experimental)
    # This is the exact transformation from cinnabar.plotting.plot_DGs
    computational_mean = dg_computational.mean()
    experimental_mean = dg_experimental.mean()
    shift = experimental_mean - computational_mean

    return dg_computational + shift


def _extract_system_edges_and_nodes(
    results: dict[str, Any],
    keys: list[str],
) -> tuple[list[dict], dict[str, list[str]], dict[str, dict], dict[str, list], dict[str, list]]:
    """Extract per-system edges and nodes from FEMaps across all keys.

    This implements Phase 2 extraction:
    - Joins computational to experimental via the cinnabar pattern
    - Fixes edge orientation per system from the first key's experimental edges
    - Computes intersection of edges/nodes across keys (D2: ragged coverage)
    - Applies ΔG centralizing shift (D3) per (system, key)
    - Records dropped edges/nodes excluded from cross-key statistics

    Parameters
    ----------
    results : dict[str, Any]
        Mapping of label → BenchmarkResults (from Phase 1).
    keys : list[str]
        Plot labels in config order.

    Returns
    -------
    systems : list[dict]
        System records: [{"system_group": str, "system_name": str, "id": "<group>/<name>"}, ...]
    missing : dict[str, list[str]]
        Keys not covering a system at all: {system_id: [key, ...]}
    dropped : dict[str, dict]
        Edges/nodes excluded from intersection: {system_id: {"edges": [edge_id, ...], "nodes": [ligand, ...]}}
    edges : dict[str, list]
        Per-system edge records with per_key data: {system_id: [{"edge_id": ..., "ligand_a": ..., ...}, ...]}
    nodes : dict[str, list]
        Per-system node records with centralized per_key DG: {system_id: [{"ligand": ..., "dg_exp": ..., ...}, ...]}
    """
    # Collect all systems across all keys
    all_systems = set()
    for result in results.values():
        all_systems.update(result.ddg_femaps().keys())

    systems = sorted(
        [{"system_group": sg, "system_name": sn, "id": f"{sg}/{sn}"} for sg, sn in all_systems],
        key=lambda x: x["id"]
    )

    missing = defaultdict(list)
    dropped = {}
    edges = {}
    nodes = {}

    # Fixed edge orientations per system (set by first key covering it)
    system_edge_orientations = {}

    for system_dict in systems:
        system_group = system_dict["system_group"]
        system_name = system_dict["system_name"]
        system_id = system_dict["id"]

        # Collect edges and nodes per key for this system
        system_edges_by_key = {}
        system_nodes_by_key = {}

        for key in keys:
            ddg_maps = results[key].ddg_femaps(source=key)
            dg_maps = results[key].dg_femaps(source=key, force_source_update=True)

            if (system_group, system_name) not in ddg_maps:
                missing[system_id].append(key)
                continue

            femap_ddg = ddg_maps[(system_group, system_name)]
            femap_dg = dg_maps[(system_group, system_name)]

            # Extract edges (relative/DDG) using cinnabar's join pattern
            df_rel = femap_ddg.get_relative_dataframe()

            # Split computational and experimental
            df_comp = df_rel[df_rel["computational"] == True].copy()
            df_exp = df_rel[df_rel["computational"] == False].copy()

            if df_exp.empty:
                logger.warning(f"{system_id} / {key}: no experimental DDG data, skipping")
                missing[system_id].append(key)
                continue

            # Join on (labelA, labelB)
            df_comp_indexed = df_comp.set_index(["labelA", "labelB"])
            df_exp_indexed = df_exp.set_index(["labelA", "labelB"])

            # Inner join: keep only edges with both computational and experimental
            joined = df_comp_indexed.join(df_exp_indexed, how="inner", lsuffix="_comp", rsuffix="_exp")

            if joined.empty:
                logger.warning(f"{system_id} / {key}: no overlapping DDG edges after join, skipping")
                missing[system_id].append(key)
                continue

            # Fix edge orientation for this system (first key sets it, others must match)
            edge_set = set(joined.index)

            if system_id not in system_edge_orientations:
                # First key covering this system: set orientation from experimental
                system_edge_orientations[system_id] = edge_set
                logger.debug(f"{system_id}: fixed orientation from {key} ({len(edge_set)} edges)")
            else:
                # Subsequent key: assert orientation matches
                expected_orientation = system_edge_orientations[system_id]
                if edge_set != expected_orientation:
                    raise ValueError(
                        f"{system_id} / {key}: edge orientation mismatch. "
                        f"Expected {len(expected_orientation)} edges matching first key's orientation, "
                        f"got {len(edge_set)} edges. This violates the plan's orientation precondition. "
                        f"Symmetric difference: {edge_set.symmetric_difference(expected_orientation)}"
                    )

            # Store edges for this key
            edges_this_key = {}
            for (ligand_a, ligand_b), row in joined.iterrows():
                edge_id = f"{ligand_a}|{ligand_b}"
                edges_this_key[edge_id] = {
                    "ligand_a": ligand_a,
                    "ligand_b": ligand_b,
                    "ddg": float(row["DDG (kcal/mol)_comp"]),
                    "unc": float(row["uncertainty (kcal/mol)_comp"]) if pd.notna(row["uncertainty (kcal/mol)_comp"]) else None,
                    "ddg_exp": float(row["DDG (kcal/mol)_exp"]),
                    "ddg_exp_unc": float(row["uncertainty (kcal/mol)_exp"]) if pd.notna(row["uncertainty (kcal/mol)_exp"]) else None,
                }

            system_edges_by_key[key] = edges_this_key

            # Extract nodes (absolute/DG)
            df_abs = femap_dg.get_absolute_dataframe()

            df_abs_comp = df_abs[df_abs["computational"] == True].copy()
            df_abs_exp = df_abs[df_abs["computational"] == False].copy()

            if df_abs_exp.empty:
                logger.warning(f"{system_id} / {key}: no experimental DG data, skipping nodes")
                continue

            df_abs_comp_indexed = df_abs_comp.set_index("label")
            df_abs_exp_indexed = df_abs_exp.set_index("label")

            joined_abs = df_abs_comp_indexed.join(df_abs_exp_indexed, how="inner", lsuffix="_comp", rsuffix="_exp")

            if joined_abs.empty:
                logger.warning(f"{system_id} / {key}: no overlapping DG nodes after join")
                continue

            # Apply centralizing shift (D3)
            dg_comp = joined_abs["DG (kcal/mol)_comp"]
            dg_exp = joined_abs["DG (kcal/mol)_exp"]
            dg_centralized = _apply_dg_centralizing_shift(dg_comp, dg_exp)

            # Store nodes for this key
            nodes_this_key = {}
            for ligand, row in joined_abs.iterrows():
                nodes_this_key[ligand] = {
                    "dg": float(dg_centralized.loc[ligand]),
                    "unc": float(row["uncertainty (kcal/mol)_comp"]) if pd.notna(row["uncertainty (kcal/mol)_comp"]) else None,
                    "dg_exp": float(row["DG (kcal/mol)_exp"]),
                    "dg_exp_unc": float(row["uncertainty (kcal/mol)_exp"]) if pd.notna(row["uncertainty (kcal/mol)_exp"]) else None,
                }

            system_nodes_by_key[key] = nodes_this_key

        # Compute intersection of edges/nodes across keys (D2: ragged coverage)
        keys_covering_system = [k for k in keys if k not in missing.get(system_id, [])]

        if not keys_covering_system:
            # No keys cover this system at all (already in missing)
            continue

        # Edge intersection
        edge_sets = [set(system_edges_by_key[k].keys()) for k in keys_covering_system if k in system_edges_by_key]
        if edge_sets:
            edge_intersection = set.intersection(*edge_sets)
        else:
            edge_intersection = set()

        # Node intersection
        node_sets = [set(system_nodes_by_key[k].keys()) for k in keys_covering_system if k in system_nodes_by_key]
        if node_sets:
            node_intersection = set.intersection(*node_sets)
        else:
            node_intersection = set()

        # Record dropped edges/nodes
        all_edges_union = set.union(*edge_sets) if edge_sets else set()
        all_nodes_union = set.union(*node_sets) if node_sets else set()

        dropped_edges = sorted(all_edges_union - edge_intersection)
        dropped_nodes = sorted(all_nodes_union - node_intersection)

        if dropped_edges or dropped_nodes:
            dropped[system_id] = {
                "edges": dropped_edges,
                "nodes": dropped_nodes,
            }

        # Build final edge records (keep each key's full data, not just intersection)
        edges_for_system = []
        all_edge_ids = sorted(all_edges_union)

        for edge_id in all_edge_ids:
            ligand_a, ligand_b = edge_id.split("|")

            # Get experimental values from first key that has this edge
            ddg_exp = None
            ddg_exp_unc = None
            for k in keys_covering_system:
                if k in system_edges_by_key and edge_id in system_edges_by_key[k]:
                    ddg_exp = system_edges_by_key[k][edge_id]["ddg_exp"]
                    ddg_exp_unc = system_edges_by_key[k][edge_id]["ddg_exp_unc"]
                    break

            per_key = {}
            for k in keys_covering_system:
                if k in system_edges_by_key and edge_id in system_edges_by_key[k]:
                    edge_data = system_edges_by_key[k][edge_id]
                    per_key[k] = {
                        "ddg": edge_data["ddg"],
                        "unc": edge_data["unc"],
                    }

            edges_for_system.append({
                "edge_id": edge_id,
                "ligand_a": ligand_a,
                "ligand_b": ligand_b,
                "ddg_exp": ddg_exp,
                "ddg_exp_unc": ddg_exp_unc,
                "per_key": per_key,
            })

        edges[system_id] = edges_for_system

        # Build final node records
        nodes_for_system = []
        all_node_ids = sorted(all_nodes_union)

        for ligand in all_node_ids:
            # Get experimental values from first key that has this node
            dg_exp = None
            dg_exp_unc = None
            for k in keys_covering_system:
                if k in system_nodes_by_key and ligand in system_nodes_by_key[k]:
                    dg_exp = system_nodes_by_key[k][ligand]["dg_exp"]
                    dg_exp_unc = system_nodes_by_key[k][ligand]["dg_exp_unc"]
                    break

            per_key = {}
            for k in keys_covering_system:
                if k in system_nodes_by_key and ligand in system_nodes_by_key[k]:
                    node_data = system_nodes_by_key[k][ligand]
                    per_key[k] = {
                        "dg": node_data["dg"],  # centralized
                        "unc": node_data["unc"],
                    }

            nodes_for_system.append({
                "ligand": ligand,
                "dg_exp": dg_exp,
                "dg_exp_unc": dg_exp_unc,
                "per_key": per_key,
            })

        nodes[system_id] = nodes_for_system

    return systems, dict(missing), dropped, edges, nodes


def _filter_femap_to_edges(femap, edge_ids: set[str]) -> Any:
    """Filter a relative FEMap to only include edges in the given set.

    Preserves all experimental measurements (computational=False) since
    cinnabar requires experimental values for statistics computation.

    Parameters
    ----------
    femap : FEMap
        The FEMap to filter.
    edge_ids : set[str]
        Set of edge IDs in "ligand_a|ligand_b" format.

    Returns
    -------
    FEMap
        A new FEMap containing only measurements for the given edges,
        plus all experimental measurements.
    """

    filtered_femap = FEMap()

    # Iterate through all measurements in the FEMap
    for measurement in femap:
        # Always keep experimental measurements
        if not measurement.computational:
            filtered_femap.add_measurement(measurement)
            continue

        # For computational measurements, check if edge matches (in either orientation)
        edge_id_forward = f"{measurement.labelA}|{measurement.labelB}"
        edge_id_reverse = f"{measurement.labelB}|{measurement.labelA}"

        if edge_id_forward in edge_ids or edge_id_reverse in edge_ids:
            filtered_femap.add_measurement(measurement)

    return filtered_femap


def _filter_femap_to_nodes(femap, node_ids: set[str]) -> Any:
    """Filter an absolute FEMap to only include nodes in the given set.

    Preserves all experimental measurements (computational=False) since
    cinnabar requires experimental values for statistics computation.

    Parameters
    ----------
    femap : FEMap
        The FEMap to filter.
    node_ids : set[str]
        Set of ligand names.

    Returns
    -------
    FEMap
        A new FEMap containing only measurements for the given nodes,
        plus all experimental measurements.
    """

    filtered_femap = FEMap()

    # Iterate through all measurements in the FEMap
    for measurement in femap:
        # Always keep experimental measurements
        if not measurement.computational:
            filtered_femap.add_measurement(measurement)
            continue

        # For absolute computational measurements, one label is ReferenceState, the other is the ligand
        # Check both labelA and labelB in case of different orientations
        ligand_label = None
        if isinstance(measurement.labelA, ReferenceState):
            ligand_label = measurement.labelB
        elif isinstance(measurement.labelB, ReferenceState):
            ligand_label = measurement.labelA
        else:
            # Both are ligands - this is actually a relative measurement
            # Include it if either ligand is in our set
            if measurement.labelA in node_ids or measurement.labelB in node_ids:
                filtered_femap.add_measurement(measurement)
            continue

        if ligand_label in node_ids:
            filtered_femap.add_measurement(measurement)

    return filtered_femap


def _compute_statistics(
    results: dict[str, Any],
    keys: list[str],
    systems: list[dict],
    edges: dict[str, list],
    nodes: dict[str, list],
    dropped: dict[str, dict],
    num_bootstraps: int,
    ecdf_bootstraps: int,
    compute_ecdf_ci: bool,
    reference_model: str | None,
) -> tuple[dict, dict]:
    """Compute cross-force-field statistics per system.

    Makes 4 compare_and_rank_results calls per system (D4: each metric
    ranked on itself). Computes ECDF with optional bootstrap CI.
    Handles dg_rmse special case (D3: centralized values, hand-rolled).

    Parameters
    ----------
    results : dict[str, Any]
        Mapping of label → BenchmarkResults.
    keys : list[str]
        Plot labels in config order.
    systems : list[dict]
        System records from Phase 2.
    edges : dict[str, list]
        Edge records from Phase 2.
    nodes : dict[str, list]
        Node records from Phase 2 (with centralized DG).
    dropped : dict[str, dict]
        Dropped edges/nodes from Phase 2.
    num_bootstraps : int
        Number of bootstrap samples for statistics.
    ecdf_bootstraps : int
        Number of bootstrap samples for ECDF CI. 0 disables CI.
    compute_ecdf_ci : bool
        Whether to compute ECDF confidence intervals.
    reference_model : str | None
        Model label used to mark bars significantly different from a reference.

    Returns
    -------
    bars : dict
        {system_id: {metric: [{"key": ..., "value": ..., "ci_lower": ..., "ci_upper": ..., "cld": ..., "different_from_reference": bool}]}}
    ecdf : dict
        {system_id: {key: {"x": [...], "y": [...], "ci_lower": [...] or None, "ci_upper": [...] or None}}}
    """
    bars = {}
    ecdf = {}

    for system_dict in systems:
        system_group = system_dict["system_group"]
        system_name = system_dict["system_name"]
        system_id = system_dict["id"]

        logger.info(f"Computing statistics for {system_id}")

        # Get intersection edges/nodes for this system
        dropped_this_system = dropped.get(system_id, {})
        dropped_edge_ids = set(dropped_this_system.get("edges", []))
        dropped_node_ids = set(dropped_this_system.get("nodes", []))

        # Build intersection edge set and node set
        edges_this_system = edges.get(system_id, [])
        nodes_this_system = nodes.get(system_id, [])

        # Intersection edges: exclude dropped
        intersection_edges = [e for e in edges_this_system if e["edge_id"] not in dropped_edge_ids]
        intersection_nodes = [n for n in nodes_this_system if n["ligand"] not in dropped_node_ids]

        # Build combined FEMaps from intersection
        # For bars/ECDF we need the FEMaps restricted to intersection
        keys_covering = [k for k in keys if any(k in e["per_key"] for e in intersection_edges)]

        if not keys_covering:
            logger.warning(f"{system_id}: no keys cover this system, skipping statistics")
            continue

        # Build intersection edge/node ID sets for filtering
        intersection_edge_ids = {e["edge_id"] for e in intersection_edges}
        intersection_node_ids = {n["ligand"] for n in intersection_nodes}

        # Build combined FEMap for DDG (intersection only)
        femap_list_ddg = []
        for key in keys_covering:
            femap_ddg = results[key].ddg_femaps(source=key)[(system_group, system_name)]

            # Filter to intersection edges only
            filtered_femap_ddg = _filter_femap_to_edges(femap_ddg, intersection_edge_ids)
            femap_list_ddg.append(filtered_femap_ddg)

        combined_femap_ddg = sum(femap_list_ddg[1:], femap_list_ddg[0]) if len(femap_list_ddg) > 1 else femap_list_ddg[0]

        # Build combined FEMap for DG (intersection only)
        femap_list_dg = []
        for key in keys_covering:
            femap_dg = results[key].dg_femaps(source=key, force_source_update=True)[(system_group, system_name)]

            # Filter to intersection nodes only
            filtered_femap_dg = _filter_femap_to_nodes(femap_dg, intersection_node_ids)
            femap_list_dg.append(filtered_femap_dg)

        combined_femap_dg = sum(femap_list_dg[1:], femap_list_dg[0]) if len(femap_list_dg) > 1 else femap_list_dg[0]

        # Assert non-ragged (Phase 2 intersection should guarantee this)
        # compare_and_rank_results will raise if ragged, but we assert here for clarity
        ddg_df = combined_femap_ddg.get_relative_dataframe()
        ddg_comp = ddg_df[ddg_df["computational"] == True]
        for source in ddg_comp["source"].unique():
            count = len(ddg_comp[ddg_comp["source"] == source])
            if count != len(intersection_edges):
                raise ValueError(
                    f"{system_id}: ragged DDG coverage detected. "
                    f"Source '{source}' has {count} edges, expected {len(intersection_edges)} (intersection). "
                    f"This should not happen after Phase 2 intersection."
                )

        dg_df = combined_femap_dg.get_absolute_dataframe()
        dg_comp = dg_df[dg_df["computational"] == True]
        for source in dg_comp["source"].unique():
            count = len(dg_comp[dg_comp["source"] == source])
            if count != len(intersection_nodes):
                raise ValueError(
                    f"{system_id}: ragged DG coverage detected. "
                    f"Source '{source}' has {count} nodes, expected {len(intersection_nodes)} (intersection). "
                    f"This should not happen after Phase 2 intersection."
                )

        # Make 4 calls to compare_and_rank_results (D4: each metric ranked on itself)
        bars_this_system = {}

        # Call 1: DDG MUE
        summary_ddg_mue, comparison_ddg_mue = compare_and_rank_results(
            femap=combined_femap_ddg,
            prediction_type="edgewise",
            rank_metric="MUE",
            metrics_to_compute=["MUE", "RMSE"],
            num_bootstraps=num_bootstraps,
            confidence_level=0.95,
            alpha=0.05,
        )
        bars_this_system["ddg_mue"] = _extract_bar_records(
            summary_ddg_mue,
            "MUE",
            keys_covering,
            _extract_reference_significance(comparison_ddg_mue, keys_covering, reference_model),
        )

        # Call 2: DDG RMSE
        summary_ddg_rmse, comparison_ddg_rmse = compare_and_rank_results(
            femap=combined_femap_ddg,
            prediction_type="edgewise",
            rank_metric="RMSE",
            metrics_to_compute=["MUE", "RMSE"],
            num_bootstraps=num_bootstraps,
            confidence_level=0.95,
            alpha=0.05,
        )
        bars_this_system["ddg_rmse"] = _extract_bar_records(
            summary_ddg_rmse,
            "RMSE",
            keys_covering,
            _extract_reference_significance(comparison_ddg_rmse, keys_covering, reference_model),
        )

        # Call 3: DG KTAU
        summary_dg_ktau, comparison_dg_ktau = compare_and_rank_results(
            femap=combined_femap_dg,
            prediction_type="nodewise",
            rank_metric="KTAU",
            metrics_to_compute=["KTAU", "RMSE"],
            num_bootstraps=num_bootstraps,
            confidence_level=0.95,
            alpha=0.05,
        )
        bars_this_system["dg_ktau"] = _extract_bar_records(
            summary_dg_ktau,
            "KTAU",
            keys_covering,
            _extract_reference_significance(comparison_dg_ktau, keys_covering, reference_model),
        )

        # Call 4: DG RMSE (for CLD only; value/CI come from centralized D3)
        summary_dg_rmse, comparison_dg_rmse = compare_and_rank_results(
            femap=combined_femap_dg,
            prediction_type="nodewise",
            rank_metric="RMSE",
            metrics_to_compute=["KTAU", "RMSE"],
            num_bootstraps=num_bootstraps,
            confidence_level=0.95,
            alpha=0.05,
        )

        reference_significance_dg_rmse = _extract_reference_significance(
            comparison_dg_rmse,
            keys_covering,
            reference_model,
        )

        # Compute dg_rmse from centralized values (D3)
        dg_centralized_by_key = {}
        dg_exp_array = None

        for node in intersection_nodes:
            ligand = node["ligand"]
            if dg_exp_array is None:
                # Build experimental array (same for all keys)
                dg_exp_array = np.array([n["dg_exp"] for n in intersection_nodes])

            for key in keys_covering:
                if key in node["per_key"]:
                    if key not in dg_centralized_by_key:
                        dg_centralized_by_key[key] = []
                    dg_centralized_by_key[key].append(node["per_key"][key]["dg"])

        # Convert to arrays
        for key in dg_centralized_by_key:
            dg_centralized_by_key[key] = np.array(dg_centralized_by_key[key])

        # Compute RMSE with bootstrap
        dg_rmse_results = _compute_dg_rmse_with_bootstrap(
            dg_centralized_by_key,
            dg_exp_array,
            num_bootstraps,
        )

        # Build dg_rmse bar records, taking CLD from summary_dg_rmse
        dg_rmse_bars = []
        for key in keys_covering:
            # Get CLD from cinnabar's call
            cld_row = summary_dg_rmse[summary_dg_rmse["Model"] == key]
            cld = cld_row["CLD"].values[0] if len(cld_row) > 0 else ""

            dg_rmse_bars.append({
                "key": key,
                "value": dg_rmse_results[key]["value"],
                "ci_lower": dg_rmse_results[key]["ci_lower"],
                "ci_upper": dg_rmse_results[key]["ci_upper"],
                "cld": cld,
                "different_from_reference": reference_significance_dg_rmse.get(key, False),
            })

        bars_this_system["dg_rmse"] = dg_rmse_bars

        bars[system_id] = bars_this_system

        # Compute ECDF for each key
        ecdf_this_system = {}

        for key in keys_covering:
            # Get DDG errors for this key from intersection edges
            errors = []
            for edge in intersection_edges:
                if key in edge["per_key"]:
                    ddg_calc = edge["per_key"][key]["ddg"]
                    ddg_exp = edge["ddg_exp"]
                    errors.append(abs(ddg_calc - ddg_exp))

            errors = np.array(errors)

            # Compute ECDF with optional CI
            num_ecdf_bootstraps = ecdf_bootstraps if compute_ecdf_ci else 0
            x, y, ci_lower, ci_upper = _compute_ecdf_with_bootstrap(errors, num_ecdf_bootstraps)

            ecdf_this_system[key] = {
                "x": x.tolist(),
                "y": y.tolist(),
                "ci_lower": ci_lower.tolist() if ci_lower is not None else None,
                "ci_upper": ci_upper.tolist() if ci_upper is not None else None,
            }

        ecdf[system_id] = ecdf_this_system

    return bars, ecdf


def _extract_bar_records(
    summary_df: pd.DataFrame,
    metric: str,
    keys: list[str],
    reference_significance: dict[str, bool] | None = None,
) -> list[dict]:
    """Extract bar records from cinnabar summary DataFrame.

    Joins on Model column by value (not by index) to avoid mislabeling
    when summary_df is sorted alphabetically.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Output from compare_and_rank_results.
    metric : str
        Metric name (MUE, RMSE, KTAU).
    keys : list[str]
        Keys in config order.
    reference_significance : dict[str, bool] | None
        Mapping of model label to whether it differs significantly from the
        selected reference model.

    Returns
    -------
    list[dict]
        [{"key": str, "value": float, "ci_lower": float, "ci_upper": float, "cld": str, "different_from_reference": bool}]
    """
    records = []
    reference_significance = reference_significance or {}

    for key in keys:
        row = summary_df[summary_df["Model"] == key]

        if len(row) == 0:
            logger.warning(f"Key '{key}' not found in summary_df Model column")
            continue

        row = row.iloc[0]

        records.append({
            "key": key,
            "value": float(row[metric]),
            "ci_lower": float(row[f"{metric}_CI_Lower"]),
            "ci_upper": float(row[f"{metric}_CI_Upper"]),
            "cld": str(row["CLD"]),
            "different_from_reference": bool(reference_significance.get(key, False)),
        })

    return records


def _extract_reference_significance(
    comparison_df: pd.DataFrame,
    keys: list[str],
    reference_model: str | None,
) -> dict[str, bool]:
    """Map each model to whether it differs significantly from a reference.

    Parameters
    ----------
    comparison_df : pd.DataFrame
        Pairwise comparison results from ``compare_and_rank_results``.
    keys : list[str]
        Keys in config order.
    reference_model : str | None
        Model label used as the reference. ``None`` disables hatching.

    Returns
    -------
    dict[str, bool]
        Mapping of model label to whether it differs significantly from the
        selected reference.
    """
    significance = {key: False for key in keys}

    if reference_model is None:
      return significance

    if reference_model not in keys:
        raise ValueError(
            f"Reference model '{reference_model}' not found in configured keys: {keys}"
        )

    for key in keys:
        if key == reference_model:
            continue

        mask = (
            ((comparison_df["Model 1"] == reference_model) & (comparison_df["Model 2"] == key)) |
            ((comparison_df["Model 1"] == key) & (comparison_df["Model 2"] == reference_model))
        )

        if not mask.any():
            logger.warning(
                "No pairwise comparison found between reference '%s' and model '%s'",
                reference_model,
                key,
            )
            continue

        significance[key] = bool(comparison_df.loc[mask, "significant"].iloc[0])

    return significance


def _minify_svg(svg_string: str) -> str:
    """Minify an SVG string while preserving legibility.

    Minification steps (in order of safety):
      1. Remove XML prolog, comments, metadata/desc elements
      2. Round path/coordinate numbers to 1 decimal
      3. Collapse inter-element whitespace

    Parameters
    ----------
    svg_string : str
        Original SVG string.

    Returns
    -------
    str
        Minified SVG string starting with <svg.
    """
    # Step 1: Remove XML prolog, comments, metadata, desc
    # Remove XML prolog
    svg = re.sub(r'<\?xml[^>]*\?>', '', svg_string)
    # Remove comments
    svg = re.sub(r'<!--.*?-->', '', svg, flags=re.DOTALL)
    # Remove metadata and desc elements
    svg = re.sub(r'<metadata[^>]*>.*?</metadata>', '', svg, flags=re.DOTALL)
    svg = re.sub(r'<desc[^>]*>.*?</desc>', '', svg, flags=re.DOTALL)

    # Step 2: Round numbers to 1 decimal
    # Match floating point numbers (but not in URLs or IDs)
    def round_number(match):
        num = float(match.group(0))
        return f"{num:.1f}"

    # Round numbers in path d attribute and other numeric attributes
    # Match pattern: number with decimal point
    svg = re.sub(r'(?<=[dMmLlHhVvCcSsQqTtAa\s,])-?\d+\.\d+', round_number, svg)

    # Step 3: Collapse whitespace between tags
    svg = re.sub(r'>\s+<', '><', svg)
    svg = svg.strip()

    return svg


def _verify_minified_svg(svg_string: str) -> tuple[bool, int]:
    """Verify a minified SVG is still valid and contains paths.

    Parameters
    ----------
    svg_string : str
        Minified SVG string.

    Returns
    -------
    is_valid : bool
        True if SVG parses and contains path elements.
    path_count : int
        Number of path elements found.
    """

    try:
        root = ET.fromstring(svg_string)
        # Count path elements (need to handle namespaces)
        path_count = len(root.findall('.//{http://www.w3.org/2000/svg}path'))
        if path_count == 0:
            # Try without namespace
            path_count = len(root.findall('.//path'))

        return path_count > 0, path_count
    except ET.ParseError:
        return False, 0


def _render_structures(
    systems: list[dict],
    edges: dict[str, list],
) -> dict[str, dict]:
    """Render 2D structure depictions for all edges.

    Uses the industry_benchmarks_network ligand network per system.
    Renders two SVGs per edge (one per ligand) with mapping highlighting.

    Parameters
    ----------
    systems : list[dict]
        System records from Phase 2.
    edges : dict[str, list]
        Edge records from Phase 2.

    Returns
    -------
    structures : dict
        {system_id: {edge_id: {"svg_a": str, "svg_b": str, "n_mapped": int,
                                "n_unique_a": int, "n_unique_b": int}}}
    """

    structures = {}
    total_edges = 0
    rendered_edges = 0
    raw_bytes = 0
    minified_bytes = 0

    # Verification samples (one per system)
    verification_samples = {}

    for system_dict in systems:
        system_group = system_dict["system_group"]
        system_name = system_dict["system_name"]
        system_id = system_dict["id"]

        logger.info(f"Rendering structures for {system_id}")

        # Load ligand network
        benchmark_data = get_data_by_system_name(system_group, system_name)
        network_path = benchmark_data.ligand_networks.get("industry_benchmarks_network")

        if network_path is None:
            logger.warning(f"{system_id}: no industry_benchmarks_network found, skipping structures")
            continue

        network = openfe.LigandNetwork.from_json(file=str(network_path))

        # Index network edges by frozenset of ligand names
        network_edge_map = {}
        for edge in network.edges:
            ligand_pair = frozenset({edge.componentA.name, edge.componentB.name})
            network_edge_map[ligand_pair] = edge

        # Render structures for each result edge
        structures_this_system = {}
        edges_this_system = edges.get(system_id, [])

        for edge_record in edges_this_system:
            total_edges += 1
            edge_id = edge_record["edge_id"]
            ligand_a, ligand_b = edge_id.split("|")

            # Look up network edge
            ligand_pair = frozenset({ligand_a, ligand_b})
            network_edge = network_edge_map.get(ligand_pair)

            if network_edge is None:
                logger.warning(f"{system_id}: no network edge for {edge_id}, skipping structure")
                continue

            # Determine orientation
            # network_edge has componentA and componentB
            # We need to render in the order specified by edge_id (ligand_a, ligand_b)
            if network_edge.componentA.name == ligand_a:
                # Network orientation matches result orientation
                mol_a = network_edge.componentA
                mol_b = network_edge.componentB
                mapping_a_to_b = network_edge.componentA_to_componentB
                mapping_b_to_a = network_edge.componentB_to_componentA
            else:
                # Network orientation is reversed
                mol_a = network_edge.componentB
                mol_b = network_edge.componentA
                mapping_a_to_b = network_edge.componentB_to_componentA
                mapping_b_to_a = network_edge.componentA_to_componentB

            # Render SVGs using draw_one_molecule_mapping
            # CRITICAL: pass the atom-index dict, not the mapping object
            try:

                # Create SVG drawers
                # Use 400x400 canvas for legibility
                d2d_a = Draw.MolDraw2DSVG(400, 400)
                d2d_b = Draw.MolDraw2DSVG(400, 400)

                # Render svg_a (shows mol_a with highlighting)
                svg_a_raw = draw_one_molecule_mapping(
                    mapping_a_to_b,  # atom-index dict
                    mol_a.to_rdkit(),
                    mol_b.to_rdkit(),
                    d2d=d2d_a
                )

                # Render svg_b (shows mol_b with highlighting)
                svg_b_raw = draw_one_molecule_mapping(
                    mapping_b_to_a,  # atom-index dict
                    mol_b.to_rdkit(),
                    mol_a.to_rdkit(),
                    d2d=d2d_b
                )

                # Assert no </script in SVGs
                if "</script" in svg_a_raw or "</script" in svg_b_raw:
                    raise ValueError(f"{system_id}/{edge_id}: SVG contains '</script' - security risk")

                # Minify
                svg_a = _minify_svg(svg_a_raw)
                svg_b = _minify_svg(svg_b_raw)

                # Assert minified output starts with <svg
                if not svg_a.startswith("<svg") or not svg_b.startswith("<svg"):
                    raise ValueError(f"{system_id}/{edge_id}: minified SVG does not start with <svg")

                # Track sizes
                raw_bytes += len(svg_a_raw) + len(svg_b_raw)
                minified_bytes += len(svg_a) + len(svg_b)

                # Store structure
                structures_this_system[edge_id] = {
                    "svg_a": svg_a,
                    "svg_b": svg_b,
                    "n_mapped": len(mapping_a_to_b),
                    "n_unique_a": mol_a.to_rdkit().GetNumAtoms() - len(mapping_a_to_b),
                    "n_unique_b": mol_b.to_rdkit().GetNumAtoms() - len(mapping_a_to_b),
                }

                rendered_edges += 1

                # Keep first edge as verification sample
                if system_id not in verification_samples:
                    verification_samples[system_id] = svg_a

            except Exception as e:
                logger.warning(f"{system_id}/{edge_id}: failed to render structure: {e}")
                continue

        if structures_this_system:
            structures[system_id] = structures_this_system

    # Verify minification didn't corrupt SVGs
    verification_passed = 0
    for system_id, svg_sample in verification_samples.items():
        is_valid, path_count = _verify_minified_svg(svg_sample)
        if is_valid:
            verification_passed += 1
        else:
            logger.warning(f"Minifier verification failed for {system_id}: parse failed or no paths")

    # Log statistics
    raw_mb = raw_bytes / (1024 * 1024)
    minified_mb = minified_bytes / (1024 * 1024)
    reduction_pct = 100 * (1 - minified_bytes / raw_bytes) if raw_bytes > 0 else 0

    logger.info(
        f"structures: {rendered_edges}/{total_edges} edges rendered, "
        f"raw {raw_mb:.2f} MB → minified {minified_mb:.2f} MB ({reduction_pct:.1f}% reduction)"
    )
    logger.info(f"minifier sanity: {verification_passed}/{len(verification_samples)} systems parse, paths retained")

    return structures


def _get_asset_cache_dir() -> Path:
    """Get the fixed asset cache directory.

    Returns the repo-local cache path examples/outputs/.widget_assets/
    which is covered by .gitignore. This path is FIXED, not keyed to
    --output-dir, so the cache survives across different output directories.

    Returns
    -------
    Path
        Asset cache directory (created if it doesn't exist).
    """
    # Repository root is two levels up from this script
    repo_root = Path(__file__).parents[2]
    cache_dir = repo_root / "examples" / "outputs" / ".widget_assets"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def _fetch_asset(url: str, cache_dir: Path) -> bytes:
    """Fetch an asset from a URL, using the cache if available.

    Parameters
    ----------
    url : str
        URL to fetch.
    cache_dir : Path
        Directory for caching downloaded assets.

    Returns
    -------
    bytes
        Asset content.

    Raises
    ------
    urllib.error.URLError
        If the fetch fails.
    ValueError
        If content-type verification fails.
    """
    # Cache key from URL (use filename)
    filename = url.split("/")[-1]
    cache_path = cache_dir / filename

    if cache_path.exists():
        logger.info(f"Using cached asset: {filename}")
        return cache_path.read_bytes()

    logger.info(f"Fetching {url}")
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read()

            # Verify content-type matches expected MIME type
            content_type = response.headers.get('Content-Type', '').lower()
            file_ext = filename.split(".")[-1].lower()

            expected_types = {
                "js": ["application/javascript", "text/javascript"],
                "css": ["text/css"],
                "woff2": ["font/woff2", "application/font-woff2", "application/octet-stream"],  # woff2 sometimes served as octet-stream
            }

            if file_ext in expected_types:
                # Check if content-type matches any expected type (handle charset suffix)
                content_type_base = content_type.split(";")[0].strip()
                if content_type_base not in expected_types[file_ext]:
                    raise ValueError(
                        f"Content-Type mismatch for {filename}: got '{content_type}', "
                        f"expected one of {expected_types[file_ext]}"
                    )

        # Write to cache
        cache_path.write_bytes(content)
        logger.info(f"Cached asset: {filename} ({len(content)} bytes)")

        return content

    except urllib.error.URLError as e:
        logger.error(f"Failed to fetch {url}: {e}")
        raise


def _rewrite_font_urls_to_data_uris(css_text: str, cache_dir: Path, base_url: str) -> tuple[str, bool]:
    """Rewrite font URLs in KaTeX CSS to base64 data: URIs.

    Finds all url(...) references in the CSS, fetches the fonts,
    and replaces the URLs with base64-encoded data: URIs.

    Parameters
    ----------
    css_text : str
        KaTeX CSS content.
    cache_dir : Path
        Directory for caching downloaded fonts.
    base_url : str
        Base URL for resolving relative font URLs (e.g., "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/").

    Returns
    -------
    css_text : str
        CSS with font URLs rewritten as data: URIs (or original URLs if fetch failed).
    all_fonts_succeeded : bool
        True if all fonts were successfully fetched and inlined, False if any failed.
    """
    # Pattern: url(fonts/KaTeX_Main-Regular.woff2)
    pattern = r'url\(([^)]+\.woff2)\)'

    # Track font fetch results
    fonts_attempted = 0
    fonts_succeeded = 0

    def replace_url(match):
        nonlocal fonts_attempted, fonts_succeeded

        font_path = match.group(1)
        # Remove quotes if present
        font_path = font_path.strip("'\"")

        # Build absolute URL
        if font_path.startswith("http"):
            font_url = font_path
        else:
            font_url = base_url + font_path

        fonts_attempted += 1

        try:
            font_bytes = _fetch_asset(font_url, cache_dir)
            # Encode as base64
            b64 = base64.b64encode(font_bytes).decode("ascii")
            fonts_succeeded += 1
            # Return data: URI
            return f'url(data:font/woff2;base64,{b64})'
        except Exception as e:
            logger.warning(f"Failed to fetch font {font_url}: {e}, leaving URL as-is")
            return match.group(0)

    css_result = re.sub(pattern, replace_url, css_text)

    all_fonts_succeeded = (fonts_attempted > 0 and fonts_attempted == fonts_succeeded)

    return css_result, all_fonts_succeeded


def _fetch_assets() -> tuple[str, str, str, str, bool]:
    """Fetch d3, KaTeX JS, and KaTeX CSS from cdnjs.

    d3 failure is fatal. KaTeX failure is non-fatal (logs WARNING, returns empty strings).

    Returns
    -------
    d3_js : str
        d3.min.js content.
    katex_js : str
        katex.min.js content (empty string on failure).
    katex_css : str
        katex.min.css content with fonts inlined as data: URIs (empty string on failure).
    cache_dir_used : str
        Path to the cache directory (for logging).
    katex_success : bool
        True if KaTeX fetch succeeded, False otherwise.
    """
    cache_dir = _get_asset_cache_dir()

    # Pinned exact versions
    d3_version = "7.9.0"
    katex_version = "0.16.9"

    d3_url = f"https://cdnjs.cloudflare.com/ajax/libs/d3/{d3_version}/d3.min.js"
    katex_js_url = f"https://cdnjs.cloudflare.com/ajax/libs/KaTeX/{katex_version}/katex.min.js"
    katex_css_url = f"https://cdnjs.cloudflare.com/ajax/libs/KaTeX/{katex_version}/katex.min.css"
    katex_base_url = f"https://cdnjs.cloudflare.com/ajax/libs/KaTeX/{katex_version}/"

    # Fetch d3 (fatal on failure)
    try:
        d3_bytes = _fetch_asset(d3_url, cache_dir)
        d3_js = d3_bytes.decode("utf-8")

        # Sanity check: verify d3 version string is present
        if d3_version not in d3_js:
            raise ValueError(f"d3.min.js does not contain expected version string '{d3_version}'")

        logger.info(f"d3.min.js: {len(d3_js)} characters, version check passed")

    except Exception as e:
        logger.error(f"FATAL: failed to fetch d3 from {d3_url}: {e}")
        raise

    # Fetch KaTeX (non-fatal on failure)
    katex_success = False
    katex_js = ""
    katex_css = ""

    try:
        katex_js_bytes = _fetch_asset(katex_js_url, cache_dir)
        katex_js = katex_js_bytes.decode("utf-8")

        # Sanity check: non-empty
        if len(katex_js) == 0:
            raise ValueError("katex.min.js is empty")

        logger.info(f"katex.min.js: {len(katex_js)} characters")

        # Fetch CSS
        katex_css_bytes = _fetch_asset(katex_css_url, cache_dir)
        katex_css_raw = katex_css_bytes.decode("utf-8")

        # Sanity check: non-empty
        if len(katex_css_raw) == 0:
            raise ValueError("katex.min.css is empty")

        logger.info(f"katex.min.css: {len(katex_css_raw)} characters (before font inlining)")

        # Rewrite font URLs to data: URIs
        katex_css, all_fonts_succeeded = _rewrite_font_urls_to_data_uris(katex_css_raw, cache_dir, katex_base_url)
        logger.info(f"katex.min.css: {len(katex_css)} characters (after font inlining)")

        if not all_fonts_succeeded:
            logger.warning("WARNING: Some KaTeX fonts failed to fetch. Treating as full KaTeX failure.")
            katex_success = False
        else:
            katex_success = True

    except Exception as e:
        logger.warning(f"WARNING: KaTeX fetch failed: {e}. Widget will render labels as literal LaTeX.")
        katex_js = ""
        katex_css = ""
        katex_success = False

    return d3_js, katex_js, katex_css, str(cache_dir), katex_success


def _assemble_html(
    data: dict,
    template_path: Path,
    widget_js_path: Path,
    d3_js: str,
    katex_js: str,
    katex_css: str,
) -> str:
    """Assemble the self-contained HTML widget.

    Substitutes placeholders in the template with inlined assets and data.

    Parameters
    ----------
    data : dict
        Widget data (will be JSON-serialized and embedded).
    template_path : Path
        Path to comparison_widget_template.html.
    widget_js_path : Path
        Path to comparison_widget.js.
    d3_js : str
        d3.min.js content.
    katex_js : str
        katex.min.js content (empty string if unavailable).
    katex_css : str
        katex.min.css content with inlined fonts (empty string if unavailable).

    Returns
    -------
    str
        Complete HTML document.
    """
    # Read template
    template_html = template_path.read_text()

    # Read widget JS
    widget_js = widget_js_path.read_text()

    # Serialize data to JSON with < escaped to prevent </script> termination
    data_json = json.dumps(data, indent=2)
    data_json_escaped = data_json.replace("<", "\\u003c")

    # Build CSS block (empty if KaTeX unavailable)
    css_block = f"<style>\n{katex_css}\n</style>" if katex_css else ""

    # Build KaTeX JS block (empty if unavailable)
    katex_js_block = f"<script>\n{katex_js}\n</script>" if katex_js else ""

    # Substitute placeholders
    html = template_html
    html = html.replace("<!-- PLACEHOLDER_CSS -->", css_block)
    html = html.replace("<!-- PLACEHOLDER_DATA -->", data_json_escaped)
    html = html.replace("<!-- PLACEHOLDER_D3 -->", f"<script>\n{d3_js}\n</script>")
    html = html.replace("<!-- PLACEHOLDER_KATEX_CSS -->", "")  # CSS already in PLACEHOLDER_CSS
    html = html.replace("<!-- PLACEHOLDER_KATEX_JS -->", katex_js_block)
    html = html.replace("<!-- PLACEHOLDER_WIDGET_JS -->", f"<script>\n{widget_js}\n</script>")

    return html


@click.command()
@click.option(
    "--config",
    required=True,
    type=click.Path(dir_okay=False),
    help="Path to YAML config file mapping plot labels to submission IDs.",
)
@click.option(
    "--output-dir",
    default="examples/outputs",
    type=str,
    help="Directory for output files (resolved against repo root). Default: examples/outputs",
)
@click.option(
    "--num-bootstraps",
    default=2000,
    type=int,
    help="Number of bootstrap samples for statistics. Default: 2000",
)
@click.option(
    "--ecdf-bootstraps",
    default=1000,
    type=int,
    help="Number of bootstrap samples for ECDF confidence bands. Default: 1000",
)
@click.option(
    "--no-ecdf-ci",
    is_flag=True,
    help="Disable ECDF confidence interval bands.",
)
@click.option(
    "--no-structures",
    is_flag=True,
    help="Skip 2D structure rendering (faster, for testing).",
)
@click.option(
    "--data-only",
    is_flag=True,
    help="Write JSON only, skip asset fetch and HTML assembly (offline testing).",
)
@click.option(
    "--seed",
    default=42,
    type=int,
    help="Random seed for bootstrap reproducibility. Default: 42",
)
@click.option(
    "--dump-counts",
    is_flag=True,
    help="Print per-system edge counts for verification (Phase 2).",
)
@click.option(
    "--reference-model",
    default=None,
    type=str,
    help="Model label used as the reference for bar hatching significance.",
)
def export_comparison_widget(
    config: str,
    output_dir: str,
    num_bootstraps: int,
    ecdf_bootstraps: int,
    no_ecdf_ci: bool,
    no_structures: bool,
    data_only: bool,
    seed: int,
    dump_counts: bool,
    reference_model: str | None,
) -> None:
    """Export widget data and optional HTML.

    Parameters
    ----------
    config : str
        YAML config path.
    output_dir : str
        Directory to write outputs.
    num_bootstraps : int
        Bootstrap count for summary statistics.
    ecdf_bootstraps : int
        Bootstrap count for ECDF confidence intervals.
    no_ecdf_ci : bool
        If ``True``, omit ECDF confidence intervals.
    no_structures : bool
        If ``True``, skip ligand depiction rendering.
    data_only : bool
        If ``True``, write JSON only.
    seed : int
        Random seed for bootstrap reproducibility.
    dump_counts : bool
        If ``True``, print verification counts.
    reference_model : str | None
        Label used to mark bars as significantly different from a reference.

    Returns
    -------
    None
        Writes files to ``output_dir``.
    """
    # Use paths as-is
    config_path = Path(config)
    output_path = Path(output_dir)

    # Ensure output directory exists
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading config from {config_path}")

    # Load config YAML
    with open(config_path, "r") as f:
        config_data = yaml.safe_load(f)

    # Validate config structure
    if not isinstance(config_data, dict):
        raise TypeError(f"Config must be a dict, got {type(config_data)}")

    if "submissions" not in config_data:
        raise KeyError("Config must have top-level 'submissions' key")

    submissions = config_data["submissions"]
    if not isinstance(submissions, dict):
        raise TypeError(f"Config 'submissions' must be a dict, got {type(submissions)}")

    if reference_model is not None and reference_model not in submissions:
        raise ValueError(
            f"Reference model '{reference_model}' not found in config submissions: {list(submissions.keys())}"
        )

    logger.info(f"Loaded {len(submissions)} submission(s) from config")

    # Load BenchmarkResults for each key
    results = {}
    for label, submission_id in submissions.items():
        logger.info(f"Loading '{label}' (submission_id: {submission_id})")
        results[label] = get_benchmark_results(submission_id)

    logger.info(f"Loaded {len(results)} keys")

    # Prime FEMap caches with source labels
    for label, result in results.items():
        logger.info(f"Priming FEMaps for '{label}'")

        # Prime ddg_femaps with source
        ddg_maps = result.ddg_femaps(source=label)
        logger.info(f"  ΔΔG: {len(ddg_maps)} system(s) — {list(ddg_maps.keys())}")

        # Prime dg_femaps with source and force_source_update=True
        # (Required per 2026-08-28_debug_custom_source_close.md: calling
        # ddg_femaps(source=X) then dg_femaps(source=X) raises ValueError
        # without force_source_update=True)
        dg_maps = result.dg_femaps(source=label, force_source_update=True)
        logger.info(f"  ΔG:  {len(dg_maps)} system(s) — {list(dg_maps.keys())}")

    logger.info("Phase 1 complete: config loaded, FEMaps primed")

    # Phase 2: Extract edges and nodes
    logger.info("Phase 2: Extracting edges and nodes per system")

    keys = list(submissions.keys())
    systems, missing, dropped, edges, nodes = _extract_system_edges_and_nodes(results, keys)

    logger.info(f"Extracted {len(systems)} system(s)")

    # Verification output for --dump-counts
    if dump_counts:
        # Per-system edge counts
        system_counts = []
        for system_dict in systems:
            system_id = system_dict["id"]
            system_name_short = system_dict["system_name"]
            edge_count = len(edges.get(system_id, []))
            system_counts.append(f"{system_name_short} {edge_count}")

        total_edges = sum(len(edges.get(s["id"], [])) for s in systems)
        dropped_edge_count = sum(len(d.get("edges", [])) for d in dropped.values())
        dropped_node_count = sum(len(d.get("nodes", [])) for d in dropped.values())

        counts_str = " / ".join(system_counts)
        print(f"{counts_str} / total {total_edges}")
        print(f"dropped: {dropped_edge_count} edges, {dropped_node_count} nodes")
        print("orientation assertion: passed")

    logger.info("Phase 2 complete: edges and nodes extracted")

    # Phase 3: Compute statistics (bars and ECDF)
    logger.info("Phase 3: Computing statistics (bars and ECDF)")

    # Seed NumPy RNG for reproducibility
    np.random.seed(seed)
    logger.info(f"Set NumPy random seed to {seed}")

    bars, ecdf_data = _compute_statistics(
        results=results,
        keys=keys,
        systems=systems,
        edges=edges,
        nodes=nodes,
        dropped=dropped,
        num_bootstraps=num_bootstraps,
        ecdf_bootstraps=ecdf_bootstraps,
        compute_ecdf_ci=not no_ecdf_ci,
        reference_model=reference_model,
    )

    logger.info(f"Computed statistics for {len(bars)} system(s)")

    # Verification output for Phase 3
    if dump_counts and "jacs_set/bace" in bars:
        bace_bars = bars["jacs_set/bace"]
        print("\nPhase 3 verification (jacs_set/bace point estimates):")
        print("ddg_mue:", ", ".join(f"{b['value']:.4f}" for b in bace_bars["ddg_mue"]))
        print("ddg_rmse:", ", ".join(f"{b['value']:.4f}" for b in bace_bars["ddg_rmse"]))
        print("(Compare to notebook cell ff9e6f25: ddg_mue = 0.9637, 0.8530, 0.7835, 0.7998, 0.7878)")
        print("(Compare to notebook cell ff9e6f25: ddg_rmse = 1.2014, 1.0713, 0.9716, 0.9895, 1.0138)")

        # Also verify dg_rmse self-consistency
        bace_dg_rmse = bace_bars["dg_rmse"]
        bace_nodes_data = nodes["jacs_set/bace"]
        print("\ndg_rmse self-consistency check:")
        for bar in bace_dg_rmse:
            key = bar["key"]
            # Recompute from nodes
            dg_calc = []
            dg_exp_list = []
            for node in bace_nodes_data:
                if key in node["per_key"]:
                    dg_calc.append(node["per_key"][key]["dg"])
                    dg_exp_list.append(node["dg_exp"])
            rmse_recomputed = np.sqrt(np.mean((np.array(dg_calc) - np.array(dg_exp_list)) ** 2))
            print(f"  {key}: emitted={bar['value']:.6f}, recomputed={rmse_recomputed:.6f}, match={abs(bar['value'] - rmse_recomputed) < 1e-6}")

        # Verify ECDF invariants
        print("\nECDF invariants check:")
        for system_id, ecdf_system in ecdf_data.items():
            for key, ecdf_series in ecdf_system.items():
                y = np.array(ecdf_series["y"])
                is_nondecreasing = np.all(np.diff(y) >= 0)
                y_max_is_one = abs(y[-1] - 1.0) < 1e-10
                x_len = len(ecdf_series["x"])
                y_len = len(y)

                # Check length matches intersection edge count
                intersection_edge_count = len([e for e in edges[system_id] if key in e["per_key"] and e["edge_id"] not in dropped.get(system_id, {}).get("edges", [])])

                if not (is_nondecreasing and y_max_is_one and x_len == y_len == intersection_edge_count):
                    print(f"  {system_id}/{key}: FAIL - nondecreasing={is_nondecreasing}, y[-1]=={y[-1]:.6f}, len(x)={x_len}, len(y)={y_len}, intersection_edges={intersection_edge_count}")
                    break
            else:
                continue
            break
        else:
            print("  All ECDF series pass invariants")

    logger.info("Phase 3 complete: statistics computed")

    # Phase 4: Render 2D structure depictions
    if no_structures:
        logger.info("Phase 4: Skipping structure rendering (--no-structures)")
        structures = {}
    else:
        logger.info("Phase 4: Rendering 2D structure depictions")
        structures = _render_structures(
            systems=systems,
            edges=edges,
        )
        logger.info("Phase 4 complete: structures rendered")

    # Build complete data structure per Data Contract
    logger.info("Assembling final data structure")

    # Build key_submission_ids mapping
    key_submission_ids = {label: submission_id for label, submission_id in submissions.items()}

    # Build complete data dict
    data = {
        "generated": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "config_path": str(config_path),
        "seed": seed,
        "reference_model": reference_model,
        "keys": keys,
        "key_submission_ids": key_submission_ids,
        "systems": systems,
        "bars": bars,
        "missing": missing,
        "dropped": dropped,
        "ecdf": ecdf_data,
        "edges": edges,
        "nodes": nodes,
        "structures": structures,
        "units": "kcal/mol",
        "notes": {
            "dg_centralized": True,
            "mapping_source": "industry_benchmarks_network",
            "katex": False if data_only else True,  # Will be updated after asset fetch; false for --data-only since no fetch occurs
        },
    }

    # Phase 6: Asset fetching and HTML assembly
    if data_only:
        logger.info("Phase 6: Skipping asset fetch and HTML assembly (--data-only)")

        # Write JSON output (--data-only mode: no asset fetch, katex=false)
        json_output_path = output_path / "comparison_widget_data.json"
        logger.info(f"Writing JSON to {json_output_path}")

        with open(json_output_path, "w") as f:
            json.dump(data, f, indent=2)

        json_size_mb = json_output_path.stat().st_size / (1024 * 1024)
        logger.info(f"JSON written: {json_size_mb:.2f} MB")

        logger.info("Export complete (data-only mode)")
        return

    logger.info("Phase 6: Fetching assets and assembling HTML")

    # Fetch assets
    d3_js, katex_js, katex_css, cache_dir_used, katex_success = _fetch_assets()

    # Update notes.katex based on actual fetch result
    data["notes"]["katex"] = katex_success

    # Write JSON output (after KaTeX flag update to ensure JSON/HTML consistency)
    json_output_path = output_path / "comparison_widget_data.json"
    logger.info(f"Writing JSON to {json_output_path}")

    with open(json_output_path, "w") as f:
        json.dump(data, f, indent=2)

    json_size_mb = json_output_path.stat().st_size / (1024 * 1024)
    logger.info(f"JSON written: {json_size_mb:.2f} MB")

    # Get template and widget JS paths
    script_dir = Path(__file__).parent
    template_path = script_dir / "comparison_widget_template.html"
    widget_js_path = script_dir / "comparison_widget.js"

    if not template_path.exists():
        raise FileNotFoundError(f"Template not found: {template_path}")
    if not widget_js_path.exists():
        raise FileNotFoundError(f"Widget JS not found: {widget_js_path}")

    # Assemble HTML
    html = _assemble_html(
        data=data,
        template_path=template_path,
        widget_js_path=widget_js_path,
        d3_js=d3_js,
        katex_js=katex_js,
        katex_css=katex_css,
    )

    # Write HTML output
    html_output_path = output_path / "comparison_widget.html"
    logger.info(f"Writing HTML to {html_output_path}")

    with open(html_output_path, "w") as f:
        f.write(html)

    html_size_mb = html_output_path.stat().st_size / (1024 * 1024)
    logger.info(f"HTML written: {html_size_mb:.2f} MB")

    # Log asset sizes
    d3_kb = len(d3_js) / 1024
    katex_js_kb = len(katex_js) / 1024 if katex_js else 0
    katex_css_kb = len(katex_css) / 1024 if katex_css else 0
    widget_js_kb = len(widget_js_path.read_text()) / 1024

    logger.info(f"Inlined assets: d3={d3_kb:.1f} KB, katex_js={katex_js_kb:.1f} KB, katex_css={katex_css_kb:.1f} KB, widget_js={widget_js_kb:.1f} KB")
    logger.info(f"Asset cache: {cache_dir_used}")

    logger.info("Phase 6 complete: HTML assembled")
    logger.info("Export complete")


if __name__ == "__main__":
    export_comparison_widget()
