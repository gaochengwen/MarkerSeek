# MarkerSeek

MarkerSeek is a command-line toolkit for chloroplast nucleotide diversity analysis. It reads multiple annotated GenBank files, aligns the full plastomes with MAFFT, calculates site-wise and window-wise Pi values, summarises Pi for genes and intergenic spacers, labels high-polymorphism regions, and exports publication-style figures tuned for Nature-size layouts.

## Features

- Multiple GenBank inputs or whole directories as the primary workflow
- Automatic MAFFT alignment with reference-coordinate projection
- Sliding-window Pi table with configurable window and step size
- Gene and intergenic-spacer Pi summary table from the reference annotation
- Automatic LSC, SSC, IRb, and IRa inference with optional manual overrides
- Hotspot labelling by top percentage, top N, or explicit Pi threshold
- mVISTA-style pairwise similarity figure across all non-reference samples
- Figure export as both `PDF` and `PNG` with a 183 mm-wide layout and Arial text

## Installation

Clone the repository and install into a virtual environment:

```bash
git clone https://github.com/gaochengwen/MarkerSeek.git
cd MarkerSeek
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install .
```

Or install directly from GitHub without cloning:

```bash
python -m pip install git+https://github.com/gaochengwen/MarkerSeek.git
```

MarkerSeek expects `mafft` to be available in `PATH`, or you can pass `--mafft-bin /path/to/mafft`. Supported Python versions: 3.11, 3.12, and 3.13.

On HPC clusters that use environment modules, you can typically make MAFFT available with something like:

```bash
module load mafft-7.490
```

## Quick Start

Run MarkerSeek on every GenBank file in `test_data/`, using `Salvia_chinensis.gb` as the coordinate reference:

```bash
markerseek analyze test_data \
  --reference test_data/Salvia_chinensis.gb \
  --outdir output \
  --window 600 \
  --step 200
```

Without `--reference`, the **first input file (sorted alphabetically)** is used as the reference. The reference defines:

- the coordinate system reported in the TSVs and figures
- the gene/IGS annotation that drives feature naming and hotspot labels
- the LSC / IRb / SSC / IRa inference

After the run completes, you should see:

- `output/pi_plot.pdf` and `output/pi_plot.png` — the Pi diversity figure with hotspot labels
- `output/similarity_plot.pdf` and `output/similarity_plot.png` — mVISTA-style identity figure
- `output/pi_windows.tsv` — sliding-window Pi values
- `output/pi_features.tsv` — per-feature (gene / tRNA / rRNA / IGS) Pi values

Example with explicit Pi threshold and tighter peak labelling:

```bash
markerseek analyze test_data \
  --reference test_data/Salvia_chinensis.gb \
  --outdir output_threshold \
  --hotspot-mode threshold \
  --hotspot-value 0.04 \
  --label-mode peak-only \
```

## Key Parameters

- `--reference`: GenBank file used as the coordinate and annotation reference. Defaults to the first input.
- `--window`: sliding-window size in bp, default `600`.
- `--step`: sliding-window step in bp, default `200`.
- `--hotspot-mode`: `top-percent`, `top-n`, or `threshold`. Default `top-percent`.
- `--hotspot-value`: paired with `--hotspot-mode`. Top percentage (default `3` = top 3%), top count, or minimum Pi value.
- `--label-mode`: `peak-only`, `all`, or `none`. Default `peak-only` shows one label per cluster of consecutive hotspot windows.
- `--label-max`: maximum number of hotspot labels drawn on the Pi figure. Default: no limit.
- `--label-min-distance`: minimum midpoint spacing in bp between labeled hotspots. Default `0` (label every peak above the threshold).
- `--similarity-window`: window size in bp for the similarity figure. Default `100`.
- `--similarity-step`: step size in bp for the similarity figure. Default `20`.
- `--similarity-floor`: lower bound of the similarity y-axis (fraction in 0–1). Default `0.5`.
- `--no-similarity-plot`: skip generating `similarity_plot.{pdf,png}`.
- `--mafft-bin`: MAFFT executable or absolute path. Default `mafft`.

Manual region overrides use 1-based inclusive coordinates and require all four:

```bash
markerseek analyze input_dir \
  --reference input_dir/sample_A.gb \
  --lsc 1:84521 \
  --irb 84522:109913 \
  --ssc 109914:128312 \
  --ira 128313:151458
```

## Outputs

### `pi_windows.tsv`

One row per sliding window across the reference genome.

| Column | Description |
| --- | --- |
| `window_id` | Stable window identifier (`W0001`, `W0002`, …). |
| `start` | 1-based inclusive start position on the reference. |
| `end` | 1-based inclusive end position on the reference. |
| `midpoint` | 1-based midpoint of the window, used as the x-coordinate in the figure. |
| `pi` | Mean per-site nucleotide diversity (Nei & Li 1979) over the window's valid sites. |
| `valid_sites` | Number of alignment columns inside the window where ≥2 samples carry canonical A/C/G/T (gaps and ambiguity codes are excluded). |
| `region` | Plastome region containing the midpoint: `LSC`, `IRb`, `SSC`, or `IRa`. |
| `label_name` | Gene/IGS name with the largest overlap to the window; used as the hotspot label on the figure. |
| `is_hotspot` | `yes` if the window is selected as a hotspot under the current `--hotspot-mode`/`--hotspot-value`, otherwise `no`. |

### `pi_features.tsv`

One row per annotated feature on the reference: genes, tRNAs, rRNAs, and intergenic spacers (IGS). For multi-exon genes, only the individual exon parts are reported (e.g. `trnK-UUU_part1`, `trnK-UUU_part2`); the whole-gene span is omitted to avoid double-counting.

| Column | Description |
| --- | --- |
| `feature_id` | Unique identifier (gene name, `gene_partN` for multi-exon parts, or `IGS_NNN` for intergenic spacers). |
| `feature_type` | `gene`, `tRNA`, `rRNA`, or `igs`. |
| `parent_gene` | Parent gene symbol; for IGS rows this is `<left>|<right>` flanking gene names. |
| `label_name` | Display name used for figures and grouping (gene symbol or `<left>-<right>` for IGS). |
| `start` | 1-based inclusive start position on the reference. |
| `end` | 1-based inclusive end position on the reference. |
| `strand` | Coding strand: `1`, `-1`, or `0` (IGS / unstranded). |
| `length_bp` | Feature length in bp on the reference (handles origin-spanning features). |
| `region` | Plastome region containing the feature midpoint: `LSC`, `IRb`, `SSC`, or `IRa`. |
| `pi` | Mean per-site nucleotide diversity (Nei & Li 1979) over the feature's valid sites; `NA` if no valid sites. |

### Figures

- `pi_plot.{pdf,png}` — sliding-window Pi curve with LSC/IRb/SSC/IRa colour bands, gene-strand track, and hotspot peak labels.
- `similarity_plot.{pdf,png}` — mVISTA-style pairwise identity tracks for every non-reference sample against the reference, with shared LSC/IRb/SSC/IRa colour bands and gene track.

## Development

```bash
python3 -m pip install -e ".[dev]"
pytest
```
