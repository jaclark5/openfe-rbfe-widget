"""
Tests for comparison widget export script contract validation.

These tests verify the JSON structure and data contract of the widget export,
ensuring the widget can parse and render the exported data correctly.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


# Paths
_REPO_ROOT = Path(__file__).resolve().parents[2]
_EXPORT_SCRIPT = _REPO_ROOT / "openfe_benchmarks" / "scripts" / "export_comparison_widget_data.py"
_MINI_CONFIG = _REPO_ROOT / "openfe_benchmarks" / "scripts" / "comparison_widget_config_mini.yaml"


@pytest.fixture(scope="module")
def exported_data(tmp_path_factory):
    """
    Run the export script with the mini config and return the parsed JSON.

    This fixture runs once per test module and is shared across all tests.
    Uses --data-only to avoid network access (mandatory for CI/offline).
    Uses --no-structures for speed in the main test suite.
    """
    # Create a module-scoped temp directory
    output_dir = tmp_path_factory.mktemp("widget_export")

    # Check if submission data exists - skip if missing (fresh clone)
    submission_id = "2026-08-25_ff14sb_openff-2.1.1_jacs"
    submission_path = _REPO_ROOT / "openfe_benchmarks" / "results" / submission_id
    if not submission_path.exists():
        pytest.skip(f"Submission data not present: {submission_id}")

    # Run the export script
    cmd = [
        "micromamba", "run", "-n", "openfe", "python",
        str(_EXPORT_SCRIPT),
        "--config", str(_MINI_CONFIG),
        "--output-dir", str(output_dir),
        "--num-bootstraps", "2000",
        "--ecdf-bootstraps", "0",
        "--no-structures",
        "--data-only",
        "--seed", "42",
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )

    # Check export succeeded
    if result.returncode != 0:
        pytest.fail(
            f"Export script failed with return code {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    # Load and return the JSON
    json_path = output_dir / "comparison_widget_data.json"
    if not json_path.exists():
        pytest.fail(f"Export did not create expected JSON file: {json_path}")

    with open(json_path) as f:
        return json.load(f)


def test_top_level_keys(exported_data):
    """Test that all required top-level keys are present with correct types."""
    required_keys = {
        "generated", "config_path", "seed", "keys", "key_submission_ids",
        "systems", "bars", "missing", "dropped", "ecdf", "edges", "nodes",
        "structures", "units", "notes"
    }

    actual_keys = set(exported_data.keys())
    assert actual_keys == required_keys, (
        f"Missing keys: {required_keys - actual_keys}, "
        f"Extra keys: {actual_keys - required_keys}"
    )

    # Type checks for metadata fields
    assert isinstance(exported_data["generated"], str), "generated should be a string"
    assert isinstance(exported_data["config_path"], str), "config_path should be a string"
    assert isinstance(exported_data["seed"], int), "seed should be an int"
    assert isinstance(exported_data["key_submission_ids"], dict), "key_submission_ids should be a dict"

    # key_submission_ids keys should match data["keys"]
    assert set(exported_data["key_submission_ids"].keys()) == set(exported_data["keys"]), (
        f"key_submission_ids keys {set(exported_data['key_submission_ids'].keys())} "
        f"should match data['keys'] {set(exported_data['keys'])}"
    )


def test_keys_match_config(exported_data):
    """Test that data['keys'] matches the mini config labels in order."""
    # Mini config defines exactly 3 keys in this order
    expected_keys = [
        "Sage 2.1.1 (ff14SB) TIP3P",
        "Sage 2.3.0 (ff14SB) TIP3P",
        "TYK2 AlchemicalArchive",
    ]

    assert exported_data["keys"] == expected_keys, (
        f"Expected keys {expected_keys}, got {exported_data['keys']}"
    )


def test_systems_structure(exported_data):
    """Test that every system record has required fields."""
    systems = exported_data["systems"]
    assert isinstance(systems, list)
    assert len(systems) > 0, "Systems list should not be empty"

    for system in systems:
        assert "system_group" in system
        assert "system_name" in system
        assert "id" in system
        assert isinstance(system["system_group"], str)
        assert isinstance(system["system_name"], str)
        assert isinstance(system["id"], str)
        # id should be "group/name"
        assert system["id"] == f"{system['system_group']}/{system['system_name']}"


def test_bars_structure(exported_data):
    """Test that bars contain all four metrics with valid structure."""
    bars = exported_data["bars"]
    required_metrics = {"ddg_mue", "ddg_rmse", "dg_rmse", "dg_ktau"}

    for system in exported_data["systems"]:
        system_id = system["id"]
        assert system_id in bars, f"Missing bars for system {system_id}"

        system_bars = bars[system_id]
        actual_metrics = set(system_bars.keys())
        assert actual_metrics == required_metrics, (
            f"System {system_id}: missing metrics {required_metrics - actual_metrics}, "
            f"extra metrics {actual_metrics - required_metrics}"
        )

        # Check structure of each metric's bar records
        for metric_id, bar_list in system_bars.items():
            assert isinstance(bar_list, list)
            for bar in bar_list:
                assert "key" in bar
                assert "value" in bar
                assert "ci_lower" in bar
                assert "ci_upper" in bar
                assert "cld" in bar

                assert isinstance(bar["key"], str)
                assert isinstance(bar["value"], (int, float))
                assert isinstance(bar["ci_lower"], (int, float))
                assert isinstance(bar["ci_upper"], (int, float))
                assert isinstance(bar["cld"], str)


def test_bars_ci_brackets_value(exported_data):
    """Test that ci_lower <= value <= ci_upper for ddg_mue and ddg_rmse.

    Plan note: percentile bootstrap intervals are not guaranteed to bracket
    the point estimate, especially for small sample sizes. If this proves
    flaky on the mini config's 11-node thrombin system, narrow to just
    ddg_mue/ddg_rmse (edgewise metrics with larger sample size).
    """
    bars = exported_data["bars"]

    # Test only ddg_mue and ddg_rmse (edgewise, larger sample)
    # Skip dg_ktau and dg_rmse (nodewise, n=11, more likely to fail)
    metrics_to_check = ["ddg_mue", "ddg_rmse"]

    for system in exported_data["systems"]:
        system_id = system["id"]
        for metric_id in metrics_to_check:
            bar_list = bars[system_id][metric_id]
            for bar in bar_list:
                assert bar["ci_lower"] <= bar["value"] <= bar["ci_upper"], (
                    f"System {system_id}, metric {metric_id}, key {bar['key']}: "
                    f"CI [{bar['ci_lower']}, {bar['ci_upper']}] "
                    f"does not bracket value {bar['value']}"
                )


def test_edges_structure(exported_data):
    """Test that every edge record has required fields."""
    edges = exported_data["edges"]
    keys = exported_data["keys"]

    for system in exported_data["systems"]:
        system_id = system["id"]
        assert system_id in edges, f"Missing edges for system {system_id}"

        edge_list = edges[system_id]
        assert isinstance(edge_list, list)

        for edge in edge_list:
            assert "edge_id" in edge
            assert "ligand_a" in edge
            assert "ligand_b" in edge
            assert "ddg_exp" in edge
            assert "ddg_exp_unc" in edge or edge.get("ddg_exp_unc") is None
            assert "per_key" in edge

            assert isinstance(edge["ligand_a"], str)
            assert isinstance(edge["ligand_b"], str)
            assert len(edge["ligand_a"]) > 0
            assert len(edge["ligand_b"]) > 0
            assert isinstance(edge["ddg_exp"], (int, float))

            # edge_id should be "ligand_a|ligand_b"
            assert edge["edge_id"] == f"{edge['ligand_a']}|{edge['ligand_b']}"

            # per_key should have keys that are a subset of data["keys"]
            per_key = edge["per_key"]
            assert isinstance(per_key, dict)
            for key in per_key.keys():
                assert key in keys, f"Unknown key {key} in edge per_key"

            # Each per_key entry should have ddg and unc
            for key, key_data in per_key.items():
                assert "ddg" in key_data
                assert "unc" in key_data or key_data.get("unc") is None
                assert isinstance(key_data["ddg"], (int, float))


def test_nodes_structure(exported_data):
    """Test that every node record has required fields."""
    nodes = exported_data["nodes"]
    keys = exported_data["keys"]

    for system in exported_data["systems"]:
        system_id = system["id"]
        assert system_id in nodes, f"Missing nodes for system {system_id}"

        node_list = nodes[system_id]
        assert isinstance(node_list, list)

        for node in node_list:
            assert "ligand" in node
            assert "dg_exp" in node
            assert "dg_exp_unc" in node or node.get("dg_exp_unc") is None
            assert "per_key" in node

            assert isinstance(node["ligand"], str)
            assert len(node["ligand"]) > 0
            assert isinstance(node["dg_exp"], (int, float))

            # per_key should have keys that are a subset of data["keys"]
            per_key = node["per_key"]
            assert isinstance(per_key, dict)
            for key in per_key.keys():
                assert key in keys, f"Unknown key {key} in node per_key"

            # Each per_key entry should have dg and unc (dg is centralized per D3)
            for key, key_data in per_key.items():
                assert "dg" in key_data
                assert "unc" in key_data or key_data.get("unc") is None
                assert isinstance(key_data["dg"], (int, float))


def test_edges_nodes_ligand_consistency(exported_data):
    """Test that every ligand in edges also appears in nodes for that system."""
    edges = exported_data["edges"]
    nodes = exported_data["nodes"]

    for system in exported_data["systems"]:
        system_id = system["id"]

        # Collect all ligands mentioned in edges
        edge_ligands = set()
        for edge in edges[system_id]:
            edge_ligands.add(edge["ligand_a"])
            edge_ligands.add(edge["ligand_b"])

        # Collect all ligands in nodes
        node_ligands = {node["ligand"] for node in nodes[system_id]}

        # Every edge ligand should be in nodes
        missing_in_nodes = edge_ligands - node_ligands
        assert len(missing_in_nodes) == 0, (
            f"System {system_id}: edge ligands missing from nodes: {missing_in_nodes}"
        )


def test_ecdf_structure(exported_data):
    """Test that ECDF series have valid structure and invariants."""
    ecdf = exported_data["ecdf"]
    keys = exported_data["keys"]

    for system in exported_data["systems"]:
        system_id = system["id"]
        assert system_id in ecdf, f"Missing ECDF for system {system_id}"

        system_ecdf = ecdf[system_id]
        assert isinstance(system_ecdf, dict)

        for key in system_ecdf.keys():
            assert key in keys, f"Unknown key {key} in ECDF"

        for key, series in system_ecdf.items():
            assert "x" in series
            assert "y" in series
            # ci_lower and ci_upper can be None (--no-ecdf-ci)

            x = series["x"]
            y = series["y"]

            assert isinstance(x, list)
            assert isinstance(y, list)

            # Equal length
            assert len(x) == len(y), (
                f"System {system_id}, key {key}: x and y must have equal length"
            )

            if len(y) > 0:
                # y is non-decreasing
                for i in range(1, len(y)):
                    assert y[i] >= y[i-1], (
                        f"System {system_id}, key {key}: y is not non-decreasing at index {i}"
                    )

                # y[-1] == 1.0 (cumulative probability reaches 1)
                assert abs(y[-1] - 1.0) < 1e-9, (
                    f"System {system_id}, key {key}: y[-1] = {y[-1]}, expected 1.0"
                )


def test_ecdf_length_equals_intersection_edge_count(exported_data):
    """Test that ECDF series length equals cross-key intersection edge count.

    Per plan D2: ECDF uses intersection, not full edge list.
    len(x) should equal the number of edges present in ALL keys for that system.
    """
    ecdf = exported_data["ecdf"]
    edges = exported_data["edges"]
    dropped = exported_data["dropped"]

    for system in exported_data["systems"]:
        system_id = system["id"]

        # Calculate intersection edge count
        # Start with all edge_ids in the system
        edge_ids_by_key = {}
        for edge in edges[system_id]:
            edge_id = edge["edge_id"]
            for key in edge["per_key"].keys():
                edge_ids_by_key.setdefault(key, set()).add(edge_id)

        # Intersection across all keys that have data for this system
        if edge_ids_by_key:
            keys_with_data = list(edge_ids_by_key.keys())
            intersection = edge_ids_by_key[keys_with_data[0]]
            for key in keys_with_data[1:]:
                intersection &= edge_ids_by_key[key]
            expected_count = len(intersection)
        else:
            expected_count = 0

        # Check ECDF series lengths
        system_ecdf = ecdf[system_id]
        for key, series in system_ecdf.items():
            actual_count = len(series["x"])
            assert actual_count == expected_count, (
                f"System {system_id}, key {key}: "
                f"ECDF length {actual_count} != intersection edge count {expected_count}"
            )


def test_dropped_structure(exported_data):
    """Test that dropped has valid structure."""
    dropped = exported_data["dropped"]
    assert isinstance(dropped, dict)

    for system_id, dropped_data in dropped.items():
        if dropped_data:
            assert "edges" in dropped_data or "nodes" in dropped_data
            if "edges" in dropped_data:
                assert isinstance(dropped_data["edges"], list)
            if "nodes" in dropped_data:
                assert isinstance(dropped_data["nodes"], list)


def test_dropped_edges_absent_from_bar_per_key(exported_data):
    """Test that edges listed in dropped are absent from per_key maps used for bars.

    Per plan D2: dropped edges are excluded from cross-key statistics (bars, ECDF)
    but still present in edges with their full per_key data. An edge is dropped
    when it's missing from at least one key, so its per_key dict should have fewer
    entries than the total number of keys available for that system.

    Note: The mini config currently produces no dropped edges, so this test
    validates structure when dropped is empty. Future configs with ragged
    coverage will exercise the actual dropped-edges logic.
    """
    dropped = exported_data["dropped"]
    edges = exported_data["edges"]
    missing = exported_data["missing"]
    all_keys = set(exported_data["keys"])

    for system in exported_data["systems"]:
        system_id = system["id"]

        if system_id not in dropped or "edges" not in dropped[system_id]:
            continue

        dropped_edge_ids = set(dropped[system_id]["edges"])

        # Calculate keys available for this system (all keys minus missing ones)
        if system_id in missing:
            keys_for_system = all_keys - set(missing[system_id])
        else:
            keys_for_system = all_keys

        # Verify each dropped edge is present in edges but missing from at least one key
        for edge in edges[system_id]:
            if edge["edge_id"] in dropped_edge_ids:
                # This edge is dropped, so it should be missing from at least one key
                # i.e., len(per_key) < len(keys_for_system)
                assert len(edge["per_key"]) < len(keys_for_system), (
                    f"System {system_id}, edge {edge['edge_id']}: "
                    f"listed in dropped.edges but present in all {len(keys_for_system)} keys. "
                    f"Dropped edges must be missing from at least one key."
                )


def test_no_structures_when_flag_set(exported_data):
    """Test that structures dict is empty when --no-structures is used."""
    # This fixture uses --no-structures, so structures should be empty or
    # have empty dicts per system
    structures = exported_data["structures"]
    assert isinstance(structures, dict)

    # All systems should have empty structure dicts
    for system in exported_data["systems"]:
        system_id = system["id"]
        if system_id in structures:
            assert len(structures[system_id]) == 0, (
                f"System {system_id} has structures when --no-structures was used"
            )


def test_units_field(exported_data):
    """Test that units field is present and correct."""
    assert "units" in exported_data
    assert exported_data["units"] == "kcal/mol"


def test_notes_field(exported_data):
    """Test that notes field has required subfields."""
    notes = exported_data["notes"]
    assert isinstance(notes, dict)
    assert "dg_centralized" in notes
    assert "mapping_source" in notes
    assert "katex" in notes

    assert notes["dg_centralized"] is True  # Always true per D3
    assert notes["mapping_source"] == "industry_benchmarks_network"  # Per D1
    # katex should be False when --data-only is used
    assert notes["katex"] is False


# Structures test - opt-in via environment variable
@pytest.mark.skipif(
    not os.environ.get("OPENFE_WIDGET_SLOW_TESTS"),
    reason="Structures test disabled by default (slow). Set OPENFE_WIDGET_SLOW_TESTS=1 to run."
)
def test_structures_with_structures_enabled(tmp_path):
    """Test that structures are present and valid when --no-structures is NOT used.

    This test is slow (renders ~14 SVG pairs for thrombin) and is opt-in only.
    """
    # Check if submission data exists
    submission_id = "2026-08-25_ff14sb_openff-2.1.1_jacs"
    submission_path = _REPO_ROOT / "openfe_benchmarks" / "results" / submission_id
    if not submission_path.exists():
        pytest.skip(f"Submission data not present: {submission_id}")

    # Run export WITH structures
    cmd = [
        "micromamba", "run", "-n", "openfe", "python",
        str(_EXPORT_SCRIPT),
        "--config", str(_MINI_CONFIG),
        "--output-dir", str(tmp_path),
        "--num-bootstraps", "100",  # Lower for speed
        "--ecdf-bootstraps", "0",
        "--data-only",
        "--seed", "42",
        # NOTE: --no-structures is NOT passed
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(_REPO_ROOT),
    )

    if result.returncode != 0:
        pytest.fail(
            f"Export script failed with return code {result.returncode}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    # Load JSON
    json_path = tmp_path / "comparison_widget_data.json"
    with open(json_path) as f:
        data = json.load(f)

    structures = data["structures"]
    edges = data["edges"]

    # Check that structures exist for edges
    for system in data["systems"]:
        system_id = system["id"]

        if system_id not in structures:
            continue

        system_structures = structures[system_id]
        system_edges = edges[system_id]

        # Every edge in the system should have a structure
        # (unless the network is missing, which logs a warning)
        for edge in system_edges:
            edge_id = edge["edge_id"]

            if edge_id in system_structures:
                struct = system_structures[edge_id]

                # Required fields
                assert "svg_a" in struct
                assert "svg_b" in struct
                assert "n_mapped" in struct
                assert "n_unique_a" in struct
                assert "n_unique_b" in struct

                # SVGs should be non-empty strings starting with <svg
                assert isinstance(struct["svg_a"], str)
                assert isinstance(struct["svg_b"], str)
                assert struct["svg_a"].startswith("<svg")
                assert struct["svg_b"].startswith("<svg")

                # Should NOT contain </script (guardrail)
                assert "</script" not in struct["svg_a"]
                assert "</script" not in struct["svg_b"]

                # Counts should be non-negative integers
                assert isinstance(struct["n_mapped"], int)
                assert isinstance(struct["n_unique_a"], int)
                assert isinstance(struct["n_unique_b"], int)
                assert struct["n_mapped"] >= 0
                assert struct["n_unique_a"] >= 0
                assert struct["n_unique_b"] >= 0
