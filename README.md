# openfe-rbfe-widget

Self-contained HTML widgets for OpenFE RBFE comparison results.

- [Want to view existing results](#want-to-view-existing-results)
- [How to create your own widget](#how-to-create-your-own-widget)

## Want to view existing results

Fastest path:

1. Open any folder under `results/`.
2. Find `comparison_widget.html`.
3. Download it if needed.
4. Double-click it, or drag it into a browser.

What you can do in the widget:

- Read the 2x2 summary bar-chart panel.
- Click a system name on the x-axis to open that system.
- Click a `ΔΔG` point to highlight an edge and show the mapped ligand pair.
- Click a `ΔG` point to highlight a ligand and show one ligand depiction.
- Use pair plots to compare two force fields directly.

## How to create your own widget

Setup once from repo root:

```bash
micromamba env create -f environment.yml
```

Create a new result folder:

1. Make a folder under `results/`.
2. Add `comparison_widget_config.yaml` with:

```yaml
submissions:
  Label 1: submission_id_1
  Label 2: submission_id_2
```

3. Add a `run.sh` like this:

```bash
#!/bin/bash
REFERENCE_MODEL="Sage 2.3.0 (ff14SB) TIP3P"

micromamba run -n openfe python ../../src/export_comparison_widget_data.py \
  --config comparison_widget_config.yaml \
  --output-dir . \
  --num-bootstraps 5000 \
  --ecdf-bootstraps 2000 \
  --seed 42 \
  --reference-model "$REFERENCE_MODEL"
```

4. Run it:

```bash
cd results/<your-folder>
./run.sh
```

Files written:

- `comparison_widget.html`
- `comparison_widget_data.json`

Direct export without `run.sh`:

```bash
micromamba run -n openfe python src/export_comparison_widget_data.py \
  --config results/<your-folder>/comparison_widget_config.yaml \
  --output-dir results/<your-folder> \
  --num-bootstraps 5000 \
  --ecdf-bootstraps 2000 \
  --seed 42 \
  --reference-model "Sage 2.3.0 (ff14SB) TIP3P"
```

Rules that matter:

- `REFERENCE_MODEL` must exactly match one label in `comparison_widget_config.yaml`.
- `comparison_widget.html` is the file you share or open.
- If you change anything in `src/`, rerun `run.sh`.

Quick checks:

```bash
node --check src/comparison_widget.js
```

```bash
micromamba run -n openfe python -m pytest src/test_comparison_widget_export.py
```