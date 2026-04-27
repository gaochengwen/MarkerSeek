# MarkerSeek

MarkerSeek is a command-line toolkit for chloroplast nucleotide diversity analysis. It reads multiple annotated GenBank files, aligns the full plastomes with MAFFT, calculates site-wise and window-wise Pi values, summarizes Pi for genes and intergenic spacers, labels high-polymorphism regions, and exports a publication-style figure tuned for Nature-size layouts.

## Features

- Multiple GenBank inputs or whole directories as the primary workflow
- Automatic MAFFT alignment with reference-coordinate projection
- Sliding-window Pi table with configurable window and step size
- Gene and intergenic-spacer Pi summary table from the reference annotation
- Automatic LSC, SSC, IRb, and IRa inference with optional manual overrides
- Hotspot labelling by top percentage, top N, or explicit Pi threshold
- Figure export as both `PDF` and `PNG` with a 183 mm-wide layout and Arial text

## Installation

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install .
```

MarkerSeek expects `mafft` to be available in `PATH`, or you can pass `--mafft-bin /path/to/mafft`.
Supported Python versions: 3.11, 3.12, and 3.13.

On this cluster, you can make MAFFT available with:

```bash
module load mafft-7.490
```

## Quick Start

```bash
markerseek analyze \
  test_data/Salvia_chinensis.gb \
  test_data/Salvia_fruticosa.gb \
  test_data/Salvia_mohavensis.gb \
  --outdir output \
  --window 600 \
  --step 200
```

Example with explicit hotspot thresholding and cleaner peak-only labels:

```bash
markerseek analyze test_data \
  --outdir output_threshold \
  --window 600 \
  --step 200 \
  --hotspot-mode threshold \
  --hotspot-value 0.04 \
  --label-mode peak-only \
  --label-max 8 \
  --label-min-distance 4000
```

This creates:

- `output/pi_plot.pdf`
- `output/pi_plot.png`
- `output/pi_windows.tsv`
- `output/pi_features.tsv`

## Key Parameters

- `--reference`: choose the GenBank file used for coordinate mapping and feature annotation
- `--window`: sliding-window size in bp, default `600`
- `--step`: sliding-window step in bp, default `200`
- `--hotspot-mode`: `top-percent`, `top-n`, or `threshold`
- `--hotspot-value`: paired with `--hotspot-mode`; default is `3` for top 3%
- `--label-mode`: `peak-only`, `all`, or `none`; default is `peak-only`
- `--label-max`: maximum number of hotspot labels drawn on the figure, default is no limit
- `--label-min-distance`: minimum spacing between labeled hotspot peaks in bp
- `--mafft-bin`: MAFFT executable or absolute path

Manual region overrides use 1-based inclusive coordinates:

```bash
markerseek analyze input_dir \
  --lsc 1:84521 \
  --irb 84522:109913 \
  --ssc 109914:128312 \
  --ira 128313:151458
```

## Outputs

### `pi_windows.tsv`

Columns:

- `window_id`
- `start`
- `end`
- `midpoint`
- `pi`
- `valid_sites`
- `region`
- `label_name`
- `is_hotspot`
- `hotspot_rank`

### `pi_features.tsv`

Columns:

- `feature_id`
- `feature_type`
- `parent_gene`
- `label_name`
- `start`
- `end`
- `length_bp`
- `region`
- `pi`
- `valid_sites`
- `spans_origin`
- `is_hotspot`
- `hotspot_rank`

## Chinese Quick Notes

- 首版是命令行工具，输入多个 GenBank，输出 1 张图和 2 个表。
- 默认按参考基因组自动推断 LSC、SSC、IRb、IRa；必要时可手动传边界。
- 默认热点规则是 top 3%，图上按 gene/IGS 名称标注热点窗口。
- 后续可以直接把 `run_analysis()` 封装到 FastAPI 或在线工具。

## Development

```bash
python3 -m pip install -e ".[dev]"
pytest
```
