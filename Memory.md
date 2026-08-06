# Memory.md — PPI Network Filtering for Disease Module Detection (thesis project)

_Last updated: 2026-08-05_

## Project summary

TUM bachelor thesis project (Chair of Data Science in Systems Biomedicine,
REPO4EU project) integrating tissue-/disease-specific PPI filtering into the
`nf-core/diseasemodulediscovery` Nextflow pipeline. `/root/ba` is this
project's root (thesis repo, tracks `origin main` at
`github.com/Schansiate/bachelor_thesis`); see `CLAUDE.md` in this folder for
full theoretical background, pipeline usage, and thesis timeline. The actual
pipeline implementation lives in the separate working fork checked out at
`~/hiwi_tum/schansiate/diseasemodulediscovery`
(`github.com/Schansiate/diseasemodulediscovery`).

This memory file covers two strands of work:
1. The pipeline-side feature branches implementing tissue-specific filtering
   (in the `diseasemodulediscovery` fork).
2. Exploratory analysis notebooks/scripts kept in this thesis repo
   (`playground.ipynb`, `merge_paxdb.py`).

## 1. Pipeline branches (`diseasemodulediscovery` fork)

Work on two related feature branches that add
**tissue-specific network/seed filtering** to the pipeline:

- `tissue_specific_filtering` — main feature branch (built on top of `dev`,
  also merges in `over_representation_update` and `node_degree_distribution`
  work along the way)
- `tissue_specific_filtering-multiple_filtering_sources` — sub-branch that
  was later merged back into `tissue_specific_filtering` (commit `6fc0182`),
  extending filtering to support multiple/custom expression data sources

Commits span **2026-04-23 → 2026-07-31** (based on `dev` merge-base to tip).

## `tissue_specific_filtering-multiple_filtering_sources` branch

Focused on making the filtering step accept more than one expression-data
source and giving users a way to supply their own filtering data.

- `e1cf627` — Removed the option to initialize the pipeline with both a
  samplesheet **and** seed/network parameters at once; pipeline now only
  accepts either a samplesheet, or a seed+network pair, not a mix.
- `0c95723` / `8b6809f` — nf-core lint + pre-commit formatting cleanup.
- `54ed99c` / `a099092` — Reworked filtering statistics to be network-centered
  rather than seed-centered.
- `ed603fa` — Network channel is now only created from the `--network`
  parameter when no samplesheet is provided.
- `de2b3a8` / `3bce1bf` — Fixed/extended error handling when parsing the
  samplesheet, plus additional filtering-statistics changes.
- `8cbb479` → `2ec0abf` — Initial samplesheet integration for tissue-specific
  filtering: added filtering threshold as a pipeline parameter, added
  ID conversion via gene symbol, merged `ch_tissue_specific_network` and
  `ch_tissue_specific_seeds` into a single `ch_tissue_specific_input` channel.
- `d63e935` — First implementation of expression-value visualization.
- `8574f0f` / `f3467e4` — First implementation of allowing a **custom
  filtering file** supplied by the user, plus lint fixes.
- `21bba00` — Added custom filtering input and a new `filtering_source`
  parameter so additional filtering approaches/sources can be selected.
- `7a974b0` — Added **TCGA** as a supported filtering source.

Net diff vs. `dev` merge-base: ~537 insertions across 13 files, centered on
`bin/tissue_specific_filtering.py` (new, ~260 lines), `bin/visualize_modules.py`,
`workflows/diseasemodulediscovery.nf`, and `nextflow_schema.json`.

## `tissue_specific_filtering` branch (main feature branch)

Superset of the above branch's work (merged via `6fc0182`), plus:

- `fc8a299` — **Initial working version of tissue-specific filtering**
  (2026-04-23): new `bin/tissue_specific_filtering.py` module,
  `modules/local/tissue_specific_filtering/main.nf`, wiring into
  `workflows/diseasemodulediscovery.nf` and `subworkflows/local/gt_diamond`.
- `2bb84fc` — First implementation of ID-space mapping (for reconciling gene
  IDs between network/seeds and expression data).
- `9dc9aba` — Bug fix following ID-space mapping work.
- `2b7504d` — Added node-degree annotations to output documentation.
- `219bb57` — Moved node-degree calculation into the graphtools parser
  (`bin/graph_tool_parser.py`) rather than a separate step.
- `48800ba` / `af736f2` — Outsourced network-degree distribution plotting to
  a separate script for formatting (later commit noted as a failing YAML
  version — may need follow-up).
- `661d4c9` — Added node-degree distributions specifically for filtered
  networks.
- `f8c531e` — Updated dependency container.
- `6452890` — Merged upstream `dev`; added expression-value distribution
  plot to the MultiQC report.
- `a965bbc` — First implementation of allowing a custom filtering file
  (mirrors the sub-branch's `8574f0f`, reconciled during merges).
- `605d926` — First implementation of expression-value visualization.
- `6394810` — Merged `over_representation_update` branch in (adds
  `custom_gmt_file` param and extends over-representation analysis to run
  g:Profiler additionally on just the added module nodes — see commit
  `31b3c67`).
- `6fc0182` — **Merge of `tissue_specific_filtering-multiple_filtering_sources`**
  into this branch, bringing in TCGA support and multi-source filtering.

Net diff vs. `dev` merge-base: ~1022 insertions / 363 deletions across 37
files. Largest new/changed files:
- `bin/tissue_specific_filtering.py` (new, +277 lines) — core filtering logic
- `bin/multiqc_formatter.py` (new, +106 lines) — MultiQC report formatting
- `workflows/diseasemodulediscovery.nf` (+112/-…) — pipeline wiring
- `bin/visualize_modules.py` (+78) — visualization updates
- `bin/graph_tool_parser.py` (+48) — node-degree calc moved here
- `docs/CONTRIBUTING.md` — largely rewritten/trimmed (-202/+~)

## Open items / things to check next session

- `af736f2 "failing yaml version of network degree distribution plot"` —
  commit message implies this was left in a broken/failing state; worth
  checking if it was fixed later or still needs attention.
- A stray `.main.nf.swp` (vim swap file) was committed in both branches
  under the tissue-specific-filtering module directory — should probably be
  removed and added to `.gitignore`.
- Confirm whether `tissue_specific_filtering-multiple_filtering_sources` can
  now be deleted/archived since it's fully merged into `tissue_specific_filtering`.

## 2. Exploratory analysis in this repo (`playground.ipynb`, `merge_paxdb.py`)

`playground.ipynb` is untracked (no git history), so changes are summarized
from its current content rather than a diff. It contains 6 cells, evolving
from a toy demo into a real PaxDb/GTEx comparison analysis:

- **Cell 0** — Toy demo: loads the graph-tool `karate` example graph, assigns
  random per-vertex "expression" values, colors nodes by expression on a
  Reds colormap, and renders it with `sfdp_layout` to
  `toy_expression_graph.pdf`. Early sanity-check for the expression-coloring
  approach later added to `bin/visualize_modules.py` in the pipeline.
- **Cell 1** — First TCGA access attempt via `xenaPython`: connects to the
  UCSC Xena `gdcHub`, looks up the "GDC TCGA Lung Adenocarcinoma (LUAD)"
  cohort, and lists its `star_tpm` datasets. Precursor to the pipeline's
  `7a974b0 "added TCGA as source"` commit.
- **Cell 2** — Downloads and loads one of those TCGA TPM datasets directly
  via `requests`/`gzip` into a pandas DataFrame (genes × samples).
- **Cells 3–5** — Three iterations of the same analysis: **clustering
  tissues by cross-referencing PaxDb protein abundance against GTEx median
  gene TPM**, to check whether protein- and mRNA-level tissue expression
  profiles agree.
  - Loads `paxdb_all_tissues.tsv` (wide format, produced by
    `merge_paxdb.py`, see below) and
    `GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz`,
    joins on shared gene symbols, log2-transforms.
  - Restricts to a curated 10-tissue subset (`PAXDB_TO_GTEX` dict) with an
    unambiguous 1:1 mapping between PaxDb and GTEx tissue naming (excludes
    sub-regions like `Brain_*` or `SKIN_FIBROBLAST` that have no single
    counterpart).
  - **Cell 3**: tissue-by-tissue Spearman correlation → seaborn
    `clustermap` (diverging `vlag` colormap, source-colored row/col bars for
    PaxDb vs. GTEx) → saved to `paxdb_gtex_tissue_clustermap.pdf`.
  - **Cell 4**: same setup but with cosine similarity instead of Spearman
    correlation, `mako` colormap — a variant tried for comparison, adds NaN
    sanity checks.
  - **Cell 5**: reverts to the Spearman-correlation version from cell 3
    (essentially a clean re-run/final version), with a `print("done")` at
    the end.
  - Net takeaway: iterating on the right similarity metric and colormap to
    visualize how well protein abundance (PaxDb) and transcript level
    (GTEx) tissue profiles cluster together — evidence-gathering for using
    PaxDb as an additional/alternative filtering source in the pipeline.

`merge_paxdb.py` — standalone script (not yet part of the pipeline) that
downloads all per-tissue PaxDb integrated protein-abundance datasets for a
given organism (default `9606`, human) from `pax-db.org` (v6.1) and merges
them into one wide TSV (`id` + one column per tissue), mirroring the layout
of the GTEx `gene_median_tpm` file used by
`bin/tissue_specific_filtering.py` in the pipeline. This is the script that
produces `paxdb_all_tissues.tsv` consumed by the notebook above.

### Open items
- `playground.ipynb` is exploratory/scratch — worth deciding whether the
  PaxDb/GTEx clustering analysis (cells 3–5) should be cleaned up into a
  proper thesis analysis script, and whether cells 0–2 (toy graph, TCGA
  Xena fetch) are still needed or can be dropped once superseded by the
  pipeline's own TCGA/visualization support.
- `merge_paxdb.py` is currently a standalone helper in this repo; consider
  whether it should move into the pipeline fork's `bin/` alongside
  `tissue_specific_filtering.py` if PaxDb becomes a supported filtering
  source there.
