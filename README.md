# MarkerSeek

<img src="MarkerSeek-logo.svg" width="400">

MarkerSeek is a command-line tool for discovering plastid-genome regions that are suitable as DNA-barcoding candidates. It aligns annotated GenBank plastomes with MAFFT, estimates nucleotide diversity (Pi), identifies high-polymorphism hotspot regions, scores candidate markers using diagnostic and primer-design evidence, and writes the results as tables, figures, and JSON payloads.

In addition to the command-line implementation, MarkerSeek is also provided as a web-based tool, available at http://www.bioseqhub.cn/markerseek.

## Biological Background

Plastid genomes are widely used in plant systematics because they are usually compact, collinearly annotated, and recoverable from genome skimming data. Classical plant DNA barcoding depends on finding loci that are variable enough to separate closely related species while being flanked by conserved sequence suitable for robust PCR amplification. MarkerSeek targets that problem at whole-plastome scale.

The nucleotide diversity index Pi is the mean pairwise nucleotide difference per valid aligned site [5]. In plastid surveys, peaks in Pi mark mutational hotspots such as intergenic spacers or rapidly evolving coding intervals. A high Pi value alone is not sufficient: candidate loci should also have reliable alignment, conserved flanks, species-level resolution, low estimated misclassification risk, and a working primer pair. MarkerSeek therefore combines Pi, feature-level diagnostics, and in-silico primer evidence into a single ranked candidate-marker table.

When a dataset has no within-species replicate samples, MarkerSeek reports replicate-dependent diagnostics as `NA` rather than substituting a default value. This affects `intraspecific_divergence`, `nearest_neighbor_discrimination`, `barcoding_gap`, and `misclassification_risk`. The score formula automatically redistributes the affected weights — see "MarkerSeek Score" below.

## Installation

MarkerSeek requires Python 3.11, 3.12, or 3.13 and a working MAFFT executable [4]. MAFFT must be available on `PATH`, or supplied through `--mafft-bin`.

```bash
git clone https://github.com/gaochengwen/MarkerSeek.git
cd MarkerSeek
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install .
```

For development and tests:

```bash
python -m pip install -e ".[dev]"
pytest -q
```

Common ways to install MAFFT: `conda install -c bioconda mafft`, an HPC environment module (`module load mafft`), or an OS package supplied by the server administrator. Primer design uses `primer3-py`, which is pulled in as a Python dependency.

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

Without `--reference`, MarkerSeek uses the first sorted input file as the reference. The reference defines the coordinate system, annotation labels, gene and intergenic-spacer catalogue, and inferred LSC / IRb / SSC / IRa regions.
## Key Parameters

- `--reference`: GenBank file used as the coordinate and annotation reference. Defaults to the first input.
- `--window`: sliding-window size in bp, default `600`.
- `--step`: sliding-window step in bp, default `200`.
- `--hotspot-mode`: `top-percent`, `top-n`, or `threshold`. Default `top-percent`.
- `--hotspot-value`: paired with `--hotspot-mode`. Top percentage (default `3` = top 3%), top count, or minimum Pi value.
- `--label-mode`: `peak-only`, `all`, or `none`. Default `peak-only` shows one label per cluster of consecutive hotspot windows.
- `--label-max`: maximum number of hotspot labels drawn on the Pi figure. Default: no limit.
- `--label-min-distance`: minimum midpoint spacing in bp between labeled hotspots. Default `0` (label every peak above the threshold).
- `--similarity-window`: window size in bp for the similarity figure. Default `200`.
- `--similarity-step`: step size in bp for the similarity figure. Default `60`.
- `--similarity-floor`: lower bound of the similarity y-axis (fraction in 0–1). Default `0.5`.
- `--no-similarity-plot`: skip generating `similarity_plot.{pdf,png}`.
- `--mafft-bin`: MAFFT executable or absolute path. Default `mafft`.

Manual region boundaries can override the automatic inference:

```bash
markerseek analyze input_dir \
  --reference input_dir/sample_A.gb \
  --lsc 1:84521 \
  --irb 84522:109913 \
  --ssc 109914:128312 \
  --ira 128313:151458
```

## `markerseek analyze` Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `inputs` | required | One or more GenBank files or directories containing `.gb`, `.gbk`, or `.genbank` files. |
| `--outdir` | `output` | Output directory for TSV, figure, ZIP, and JSON payload files. |
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
| `--similarity-floor` | `0.5` | Lower y-axis bound for pairwise similarity tracks (fraction in 0–1). |
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
| `--primer-mismatch` | `1` | Allowed primer-template mismatches outside the exact 3' anchor. |
| `--primer-anchor-bp` | `5` | Number of 3' bases that must match exactly. |

## Candidate-Marker Selection

`candidate_marker_features.tsv` records every projected feature (gene, tRNA, rRNA, intergenic spacer). The per-feature `feature_payload/*.json` set is restricted to a curated subset:

* every feature whose `label_name` matches a hotspot label rendered on the Pi plot, listed in **genomic order**, plus
* the **top 10 remaining features by `markerseek_score`**.

Primer design (`--primer-design`) only runs on this same set, so each feature in the table either has a primer pair or is annotated as `primer_available = no`.

## Output Files

### `pi_windows.tsv`

| Column | Description |
| --- | --- |
| `window_id` | Stable window identifier (`W0001`, `W0002`, …). |
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
| `markerseek_score` | Composite 0–100 candidate-marker score. |

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

This 19-column table is written when `--primer-design` is enabled.

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
| `primer3_penalty` | primer3 pair penalty. |
| `target_start` | 1-based target start coordinate. |
| `target_end` | 1-based target end coordinate. |
| `amplicon_size` | In-silico amplicon length range across successfully amplified samples, written as `min-max` in base pairs. |
| `cross_species_success_rate` | Fraction of samples in which the primer pair amplifies. |
| `primer_score` | 0–100 primer-pair score. |

### Other Outputs

| File | Description |
| --- | --- |
| `primer_summary.png` | Reference-amplicon length bar chart and cross-sample amplification heatmap. |
| `pi_plot.{png,pdf}` | Sliding-window Pi plot with plastome regions and hotspot labels. |
| `similarity_plot.{png,pdf}` | mVISTA-style pairwise similarity plot against the reference. |
| `feature_payload/*.json` | Per-feature JSON payloads with local Pi curve, SNP / indel positions, haplotype network, and PCA inputs for downstream visualisation. |
| `results.zip` | Archive of the downloadable result files (PNG figures excluded — they are kept alongside the archive). |

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

Default weights:

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

`variable_site_density` and `indel_density` divide counts by `length_bp`. `flanking_conservation_min` is `min(conserved_left_bp, conserved_right_bp) / 200`. `length_suitability` is trapezoidal: 0 below 300 bp, ramps to 1 from 300–600 bp, holds at 1 from 600–1500 bp, declines to 0.5 at 3000 bp, and to 0 by 5000 bp.

### Effective weights for replicate-free datasets

`barcoding_gap` and `nearest_neighbor_discrimination` are computable only when at least one species has multiple individuals. For collections that store one reference per species — the typical NCBI plastome layout — those metrics are `NA` for every feature, and including them with weight 0.10 each would impose a structural ~20-point ceiling on the score.

When MarkerSeek detects that a metric is `NA` across the entire dataset, it drops that metric from the weight set and **re-normalizes the remaining weights so their sum stays 1.0**. The score remains a valid 0–100 quantity, scaled to the information available. The actual weight set used in a run is recorded on `result.score_weights`. Tests cover both the default-weights path and the replicate-free re-weighting path.

## Primer Design and Scoring

`--primer-design` produces primer pairs only for the curated candidate-marker set described above. For each candidate region:

1. Conserved flanking windows are identified from the genome-wide alignment.
2. primer3 designs candidate pairs against the **reference** sequence.
3. Each pair is validated by full in-silico PCR on the **reference** to confirm a unique, length-conformant amplicon.
4. Universality on the **non-reference** samples is checked with a binding-only fast path: the forward primer and the reverse-complemented reverse primer are searched (anchor-strict / body-fuzzy) and the first valid `(fwd_pos, rev_pos)` combination within `[min_amplicon, max_amplicon]` confirms amplification. The amplicon length collected from this check feeds the cross-species length statistics.
5. Per-amplicon variable / indel sites are computed by **slicing the genome-wide MAFFT alignment** at the reference amplicon coordinates — no per-pair MAFFT subprocess is spawned.

Pairs are ranked by:

```text
penalty_term = 1 - clip(primer3_penalty / 5, 0, 1)
var_density  = amplicon_variable_sites / max(amplicon_mean_len, 1)
indel_density = amplicon_indel_sites  / max(amplicon_mean_len, 1)
info_term    = 0.4 * clip(var_density  / 0.10, 0, 1)
             + 0.2 * clip(indel_density / 0.05, 0, 1)
primer_score = round(100 * (0.4 * penalty_term
                          + 0.4 * cross_species_success_rate
                          + 0.2 * info_term), 1)
```

If no sample (including the reference) yields a valid single amplicon, the primer score is 0.

## In-Silico PCR

MarkerSeek searches each ungapped sample for the forward primer and the reverse-complemented reverse primer. Matching is fuzzy outside a strict 3' anchor: by default, the terminal 5 bp at the primer 3' end must match exactly, while the rest of the primer body may contain at most one mismatch. Forward / reverse hits are paired only when the reverse hit lies downstream and the resulting amplicon length falls within the configured bounds. On the reference, all valid combinations are enumerated and the shortest is taken; multiple valid combinations mark the reference hit as ambiguous and the pair is rejected. On non-reference samples, the search short-circuits on the first valid combination — only universality and amplicon length are needed there.

## Alignment Reliability

Feature-level alignment reliability is the fraction of alignment columns that pass all three filters:

```text
gap_ratio < 0.5
ambiguity_ratio < 0.5
Shannon_entropy(ACGT frequencies) < 0.95 * log2(4)
```

This keeps high-gap, high-ambiguity, and nearly saturated columns from inflating candidate quality.

## Performance and Scaling

Genome-wide MAFFT alignment uses `--retree 2 --maxiterate 0` (FFT-NS-2, progressive only). This keeps the wall-clock close to `O(N log N)` for hundreds of plastomes and avoids the iterative-refinement path (`disttbfast`, `O(N²)`) that `mafft --auto` switches to past ~30 sequences and that can exhaust memory on small servers.

Primer design + in-silico PCR is the second-largest cost. The two structural choices that keep it tractable are: (a) primer3 evaluation only on the reference, with binding-only universality checks on the rest, short-circuiting on the first valid pair per sample; (b) per-amplicon variable / indel statistics computed by reusing the genome-wide alignment, never spawning a per-pair MAFFT subprocess.

Output writing is otherwise dominated by matplotlib figures (Pi, similarity, primer summary). The raw amplicon-FASTA outputs (`primer_amplicons.fasta`, `primer_amplicons_alignment.fasta`) were deliberately removed to drop a per-primer MAFFT subprocess loop that previously ran during `write_analysis_outputs`.

## Citation

If MarkerSeek supports your analysis, cite this repository and the underlying methods and software listed below.

## References

[1] Hebert PDN, Cywinska A, Ball SL, deWaard JR. Biological identifications through DNA barcodes. *Proceedings of the Royal Society B*. 2003;270:313-321.

[2] Meyer CP, Paulay G. DNA barcoding: error rates based on comprehensive sampling. *PLoS Biology*. 2005;3:e422.

[3] Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. Primer3: new capabilities and interfaces. *Nucleic Acids Research*. 2012;40:e115.

[4] Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*. 2013;30:772-780.

[5] Nei M, Li WH. Mathematical model for studying genetic variation in terms of restriction endonucleases. *Proceedings of the National Academy of Sciences USA*. 1979;76:5269-5273.
