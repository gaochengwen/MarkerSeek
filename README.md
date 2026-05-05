# MarkerSeek

MarkerSeek is a command-line and web application for discovering plastid genome regions that are suitable as DNA barcoding candidates. It aligns annotated GenBank plastomes, estimates nucleotide diversity, identifies high-polymorphism hotspot regions, scores candidate markers with diagnostic and practical primer-design evidence, and produces tables, figures, and interactive detail-page payloads for downstream inspection.

## Biological Background

Plastid genomes are widely used in plant systematics because they are usually compact, collinearly annotated, and recoverable from genome skimming data. Classical plant DNA barcoding depends on finding loci that are variable enough to separate closely related species while remaining flanked by conserved sequence suitable for robust PCR amplification. MarkerSeek targets that problem at whole-plastome scale.

The nucleotide diversity index, Pi, is the mean pairwise nucleotide difference per valid aligned site [5]. In plastid genome surveys, peaks in Pi often mark mutational hotspots such as intergenic spacers or rapidly evolving coding intervals. A high Pi value alone is not sufficient for a marker recommendation: candidate loci should also have reliable alignment, conserved flanking sequence, species-level resolution, low estimated misclassification risk, and practical primer pairs. MarkerSeek therefore combines Pi, feature-level diagnostics, and in-silico primer performance into a ranked candidate-marker table.

When a dataset has no within-species sample pairs, MarkerSeek reports within-species-dependent diagnostics as `NA` rather than substituting a default value. This affects `intraspecific_divergence`, `nearest_neighbor_discrimination`, `barcoding_gap`, and `misclassification_risk`.

## Installation

MarkerSeek requires Python 3.11, 3.12, or 3.13 and a working MAFFT executable [4]. MAFFT must be available on `PATH`, or supplied with `--mafft-bin`.

```bash
git clone https://github.com/gaochengwen/MarkerSeek.git
cd MarkerSeek
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install .
```

For development and testing:

```bash
python -m pip install -e ".[dev]"
pytest -q
```

Typical system dependency options include `conda install -c bioconda mafft`, an HPC environment module such as `module load mafft`, or an operating-system package supplied by the server administrator. Primer design uses `primer3-py`, which is installed as a Python dependency.

## CLI Usage

Run the bundled Salvia example data with an explicit reference:

```bash
markerseek analyze test_data \
  --reference test_data/Salvia_chinensis.gb \
  --outdir output
```

Run candidate-marker discovery with primer design:

```bash
markerseek analyze test_data \
  --reference test_data/Salvia_chinensis.gb \
  --outdir output_primer \
  --primer-design
```

Without `--reference`, MarkerSeek uses the first sorted input file as the reference. The reference defines the coordinate system, annotation labels, gene and intergenic-spacer catalogue, and inferred LSC/IRb/SSC/IRa regions.

Manual plastome region boundaries can be supplied when automatic inference is unsuitable:

```bash
markerseek analyze input_dir \
  --reference input_dir/sample_A.gb \
  --lsc 1:84521 \
  --irb 84522:109913 \
  --ssc 109914:128312 \
  --ira 128313:151458
```

## Web Usage

Start the local web server:

```bash
markerseek-web --host 0.0.0.0 --port 8000
```

Open `http://localhost:8000/markerseek`, upload GenBank files, choose a reference if needed, set analysis parameters, and submit the job. The web application returns a job ID that can be pasted into the View Results page. Successful jobs show figures, downloadable outputs, a Candidate markers table, and links to hotspot detail pages. No screenshots are bundled in this repository; the web UI is self-contained and can be inspected after starting the server.

The web server stores jobs under `markerseek_jobs/` by default and records metadata in SQLite. Useful environment variables are:

The primary navigation also includes an Example page that links to the permanent demonstration job `MSK-EXAMPLE-DEMO` when that job has been generated in the local job registry.

| Variable | Description |
| --- | --- |
| `MARKERSEEK_WEB_DATA` | Directory for uploads, job outputs, and `jobs.sqlite3`. |
| `MARKERSEEK_RETENTION_DAYS` | Retention period for ordinary jobs; default `7`. |
| `MARKERSEEK_MAX_UPLOAD_BYTES` | Total upload limit per job; default `20971520`. |
| `MARKERSEEK_MAFFT_BIN` | MAFFT executable or absolute path; default `mafft`. |

## `markerseek analyze` Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `inputs` | required | One or more GenBank files or directories containing `.gb`, `.gbk`, or `.genbank` files. |
| `--outdir` | `output` | Output directory for TSV, FASTA, figure, ZIP, and JSON payload files. |
| `--reference` | first sorted input | GenBank file used as coordinate and annotation reference. |
| `--window` | `600` | Sliding-window size in base pairs for Pi estimation. |
| `--step` | `200` | Sliding-window step size in base pairs. |
| `--hotspot-mode` | `top-percent` | Hotspot selection mode: `top-percent`, `top-n`, or `threshold`. |
| `--hotspot-value` | `3.0` | Top percentage, top count, or minimum Pi threshold, depending on `--hotspot-mode`. |
| `--mafft-bin` | `mafft` | MAFFT executable name or path. |
| `--mafft-threads` | MAFFT default | Number of MAFFT worker threads. |
| `--label-mode` | `peak-only` | Pi-plot hotspot label mode: `peak-only`, `all`, or `none`. |
| `--label-max` | unlimited | Maximum number of hotspot labels drawn on the Pi plot. |
| `--label-min-distance` | `0` | Minimum midpoint spacing in base pairs between plotted hotspot labels. |
| `--similarity-window` | `200` | Window size in base pairs for the similarity figure. |
| `--similarity-step` | `60` | Step size in base pairs for the similarity figure. |
| `--similarity-floor` | `0.5` | Lower y-axis bound for pairwise similarity tracks, as a fraction from 0 to 1. |
| `--no-similarity-plot` | disabled | Skip `similarity_plot.{png,pdf}` generation. |
| `--lsc`, `--irb`, `--ssc`, `--ira` | automatic | Manual 1-based inclusive plastome region boundaries. All four must be supplied together. |
| `--primer-design` | disabled | Design primers for hotspot candidate markers and run in-silico PCR. |
| `--primer-tm-min` | `52.0` | Minimum primer melting temperature. |
| `--primer-tm-max` | `70.0` | Maximum primer melting temperature. |
| `--primer-tm-opt` | `58.0` | Optimal primer melting temperature. |
| `--primer-len-min` | `18` | Minimum primer length. |
| `--primer-len-max` | `27` | Maximum primer length. |
| `--primer-len-opt` | `20` | Optimal primer length. |
| `--primer-amplicon-min` | `80` | Minimum allowed in-silico amplicon length. |
| `--primer-amplicon-max` | `3000` | Maximum allowed in-silico amplicon length. |
| `--primer-mismatch` | `1` | Allowed primer-template mismatches outside the exact 3-prime anchor. |
| `--primer-anchor-bp` | `5` | Number of 3-prime bases that must match exactly. |

## Output Files

### `pi_windows.tsv`

| Column | Description |
| --- | --- |
| `window_id` | Stable window identifier (`W0001`, `W0002`, ...). |
| `start` | 1-based inclusive start coordinate on the reference. |
| `end` | 1-based inclusive end coordinate on the reference. |
| `midpoint` | 1-based midpoint used for plotting. |
| `pi` | Mean nucleotide diversity across valid sites in the window. |
| `valid_sites` | Number of sites where at least two samples have canonical A/C/G/T bases. |
| `region` | Plastome region at the window midpoint: `LSC`, `IRb`, `SSC`, or `IRa`. |
| `label_name` | Reference annotation label with largest overlap to the window. |
| `is_hotspot` | `yes` when selected by the configured hotspot rule, otherwise `no`. |

### `candidate_marker_features.tsv`

This 25-column table reports genes, tRNAs, rRNAs, and intergenic spacers projected onto the reference. Multi-part genes are reported as parts to avoid double counting.

| Column | Description |
| --- | --- |
| `feature_id` | Internal stable identifier, used for detail-page routing and payload names. |
| `feature_type` | Feature class: `gene`, `tRNA`, `rRNA`, or `igs`. |
| `parent_gene` | Parent gene symbol; for IGS rows, `<left>|<right>` flanking labels. |
| `label_name` | Display label, usually a gene name or `<left>-<right>` IGS label. |
| `start` | 1-based inclusive start coordinate on the reference. |
| `end` | 1-based inclusive end coordinate on the reference. |
| `strand` | `1`, `-1`, or `0` for unstranded intervals. |
| `length_bp` | Reference-coordinate feature length in base pairs. |
| `region` | Plastome region containing the feature midpoint. |
| `pi` | Mean feature-level Pi over valid sites. |
| `variable_sites` | Count of columns with at least two canonical nucleotide states. |
| `indel_sites` | Count of mixed gap/base columns. |
| `conserved_left_bp` | Conserved contiguous bases immediately upstream of the feature. |
| `conserved_right_bp` | Conserved contiguous bases immediately downstream of the feature. |
| `primer_available` | `yes`, `no`, or `NA` depending on primer-design status and success. |
| `species_resolution` | Fraction of species whose haplotypes are not shared across species. |
| `unique_haplotype_count` | Number of distinct normalized haplotypes. |
| `species_specific_haplotype_ratio` | Fraction of haplotypes observed in only one species. |
| `interspecific_divergence` | Mean valid pairwise p-distance among different species. |
| `intraspecific_divergence` | Mean valid pairwise p-distance within species; `NA` without within-species pairs. |
| `nearest_neighbor_discrimination` | Fraction of samples whose nearest neighbor includes the same species; `NA` without within-species pairs. |
| `barcoding_gap` | Minimum interspecific distance minus maximum intraspecific distance; `NA` without within-species pairs. |
| `misclassification_risk` | `1 - nearest_neighbor_discrimination`; `NA` without within-species pairs. |
| `alignment_reliability` | Fraction of columns passing gap, ambiguity, and entropy reliability filters. |
| `markerseek_score` | Composite 0-100 candidate-marker score. |

### `haplotype_assignments.tsv`

| Column | Description |
| --- | --- |
| `feature_id` | Internal feature identifier. |
| `sample_name` | Sanitized sample name. |
| `haplotype_id` | Feature-level haplotype assignment such as `H001`. |

### `sample_metadata.tsv`

| Column | Description |
| --- | --- |
| `sample_name` | Sanitized sample name used throughout outputs. |
| `species` | Species inferred from GenBank organism metadata or filename fallback. |
| `source_path` | Input GenBank path. |

### `primers.tsv`

This 23-column table is written when `--primer-design` is enabled. The first two columns are `primer_id` and `label_name`.

| Column | Description |
| --- | --- |
| `primer_id` | Primer-pair identifier formatted as `{label_name}_p{rank}`. |
| `label_name` | Candidate-marker display label associated with the primer pair. |
| `rank` | Primer-pair rank within the candidate marker. |
| `fwd_seq` | Forward primer sequence. |
| `rev_seq` | Reverse primer sequence. |
| `fwd_len` | Forward primer length. |
| `rev_len` | Reverse primer length. |
| `fwd_gc` | Forward primer GC percentage. |
| `rev_gc` | Reverse primer GC percentage. |
| `fwd_tm` | Forward primer melting temperature. |
| `rev_tm` | Reverse primer melting temperature. |
| `fwd_self_any_th` | Forward primer self-complementarity metric from primer3. |
| `rev_self_any_th` | Reverse primer self-complementarity metric from primer3. |
| `primer3_penalty` | Primer3 pair penalty. |
| `target_start` | 1-based target start coordinate. |
| `target_end` | 1-based target end coordinate. |
| `amplicon_min_len` | Minimum successful in-silico amplicon length. |
| `amplicon_max_len` | Maximum successful in-silico amplicon length. |
| `amplicon_mean_len` | Mean successful in-silico amplicon length. |
| `cross_species_success_rate` | Fraction of samples with one valid amplicon. |
| `amplicon_variable_sites` | Variable sites in successful amplicons. |
| `amplicon_indel_sites` | Indel sites in successful amplicons. |
| `primer_score` | 0-100 primer-pair score. |

### Other Outputs

| File | Description |
| --- | --- |
| `primer_amplicons.fasta` | Successful in-silico amplicons, grouped by `primer_id` and sample. |
| `primer_amplicons_alignment.fasta` | MAFFT alignments of primer amplicon groups. |
| `primer_summary.png` | Amplicon length and cross-sample amplification summary. |
| `pi_plot.{png,pdf}` | Sliding-window Pi plot with plastome regions and hotspot labels. |
| `similarity_plot.{png,pdf}` | mVISTA-style pairwise similarity plot against the reference. |
| `feature_payload/*.json` | Internal web payloads for hotspot detail pages. |
| `results.zip` | Web-server archive containing downloadable result files. |

## MarkerSeek Score

MarkerSeek uses a normalized weighted sum:

```text
score = round(100 * Σ_i w_i * n_i, 1)
```

For higher-is-better metrics:

```text
n_i = clip((x_i - lo_i) / (hi_i - lo_i), 0, 1)
```

For lower-is-better metrics:

```text
n_i = 1 - clip((x_i - lo_i) / (hi_i - lo_i), 0, 1)
```

Missing metrics contribute `0` for their weighted term. Weights sum to 1.0.

| Metric | Direction | Low | High | Weight |
| --- | --- | ---: | ---: | ---: |
| `pi` | higher better | 0 | 0.05 | 0.18 |
| `variable_site_density` | higher better | 0 | 0.10 | 0.10 |
| `indel_density` | higher better | 0 | 0.05 | 0.05 |
| `flanking_conservation_min` | higher better | 0 | 1.0 | 0.12 |
| `missing_ambig_ratio` | lower better | 0 | 0.20 | 0.05 |
| `alignment_reliability` | higher better | 0 | 1.0 | 0.10 |
| `species_resolution` | higher better | 0 | 1.0 | 0.15 |
| `barcoding_gap` | higher better | -0.02 | 0.05 | 0.10 |
| `nearest_neighbor_discrimination` | higher better | 0 | 1.0 | 0.10 |
| `length_suitability` | higher better | 0 | 1.0 | 0.05 |

`variable_site_density` and `indel_density` divide counts by `length_bp`. `flanking_conservation_min` is `min(conserved_left_bp, conserved_right_bp) / 200`. `length_suitability` is trapezoidal: 0 below 300 bp, ramps to 1 from 300-600 bp, remains 1 from 600-1500 bp, declines to 0.5 at 3000 bp, and declines to 0 by 5000 bp.

## Primer Scoring

Primer pairs are ranked by primer3 penalty, cross-sample amplification success, and amplicon information:

```text
penalty_term = 1 - clip(primer3_penalty / 5, 0, 1)
var_density = amplicon_variable_sites / max(amplicon_mean_len, 1)
indel_density = amplicon_indel_sites / max(amplicon_mean_len, 1)
info_term = 0.4 * clip(var_density / 0.10, 0, 1)
          + 0.2 * clip(indel_density / 0.05, 0, 1)
primer_score = round(100 * (0.4 * penalty_term
                          + 0.4 * cross_species_success_rate
                          + 0.2 * info_term), 1)
```

If no sample yields a valid single amplicon, the primer score is 0.

## In-Silico PCR

MarkerSeek searches each ungapped sample sequence for the forward primer and the reverse-complemented reverse primer. Matching is fuzzy outside a strict 3-prime anchor: by default, the terminal 5 bp at the primer 3-prime end must match exactly, while the remaining primer body may contain at most one mismatch. Candidate forward/reverse hits are paired only when the reverse hit lies downstream and the resulting amplicon length is within the configured bounds. Samples with multiple valid amplicons are excluded from successful-amplification counts.

## Alignment Reliability

Feature-level alignment reliability is the fraction of alignment columns that pass all three filters:

```text
gap_ratio < 0.5
ambiguity_ratio < 0.5
Shannon_entropy(ACGT frequencies) < 0.95 * log2(4)
```

This keeps high-gap, high-ambiguity, and nearly saturated columns from inflating candidate quality.

## Hotspot Detail Pages

The web detail page for each candidate marker displays feature metadata, local Pi curve, SNP and indel positions, an alignment viewer, haplotype network, species PCA, primer ranking, and FASTA/CSV downloads. The Species PCA uses marker-only plotting with sample and species names in hover text to avoid label overlap. If within-species replicate data are absent, the relevant diagnostic fields show "Insufficient replicates" with an explanatory tooltip.

## Citation

If MarkerSeek supports your analysis, cite this repository and the underlying methods and software listed below.

## References

[1] Hebert PDN, Cywinska A, Ball SL, deWaard JR. Biological identifications through DNA barcodes. *Proceedings of the Royal Society B*. 2003;270:313-321.

[2] Meyer CP, Paulay G. DNA barcoding: error rates based on comprehensive sampling. *PLoS Biology*. 2005;3:e422.

[3] Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. Primer3: new capabilities and interfaces. *Nucleic Acids Research*. 2012;40:e115.

[4] Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*. 2013;30:772-780.

[5] Nei M, Li WH. Mathematical model for studying genetic variation in terms of restriction endonucleases. *Proceedings of the National Academy of Sciences USA*. 1979;76:5269-5273.
