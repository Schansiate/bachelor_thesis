# Memory.md — PPI Network Filtering for Disease Module Detection (thesis project)

_Last updated: 2026-09-01_

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
   (as of 2026-08-07: `pipeline_experimentation.ipynb`,
   `filtering_source_exploration.ipynb`, `filtering_results.ipynb`,
   `merge_paxdb.py` — see §2 for how these were split out).

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

## 2026-08-09/10 — per-row `source` and `threshold` samplesheet columns (`tissue_specific_filtering` branch, uncommitted)

Follow-up work on the `tissue_specific_filtering` branch (still working-tree
changes, not yet committed as of 2026-08-10). Previously `--filtering_source`
and `--filtering_threshold` were pipeline-wide params applied identically to
every samplesheet row; this made it impossible to filter different
seed/network pairs against different expression sources or thresholds in one
run. Both are now per-row samplesheet columns:

- **`assets/schema_input.json`** — added optional `source` (string) and
  `threshold` (number) columns alongside the existing `tissue` column.
- **`subworkflows/local/utils_nfcore_diseasemodulediscovery_pipeline/main.nf`**
  — samplesheet parsing now destructures 7 columns (was 5); a row that sets
  `tissue` but leaves `source` and/or `threshold` blank now hard-errors
  ("must all be specified together"), no silent fallback to the global
  params. The non-samplesheet (`--seeds`/`--network` CLI) path still uses the
  global `--filtering_source`/`--filtering_threshold` params, appended via
  `.combine(...)` to keep tuple arity consistent with the samplesheet path.
- **`workflows/diseasemodulediscovery.nf`** — `ch_tissue_specific_network`/
  `ch_tissue_specific_seeds` thread `source` and `threshold` through the
  channel tuples. Both are folded into `meta.id`/`meta.network_id` as
  `<network>.<tissue>.<source>.<threshold>` — this was necessary so that (a)
  output files across different source/threshold combinations for the same
  seed/network/tissue don't collide, and (b) the samplesheet can list the
  *same* seed/network/tissue/source combination multiple times with
  different thresholds. `network_id` is also what `ch_seeds` and
  `ch_network_gt` are joined on downstream (line ~210-211), so seeds and
  source/threshold must be kept in sync there too.
- **`modules/local/tissue_specific_filtering/main.nf`** — process input
  tuple gained `val(source)` and `val(threshold)`; the script call now uses
  `${source}`/`${threshold}` instead of `${params.filtering_source}`/
  `${params.filtering_threshold}`.
- **`bin/tissue_specific_filtering.py`** — `filter_network()` and
  `save_expression_distribution()` now take `source` and `threshold` and
  save output files as `{stem}.{tissue}.{source}.{threshold}.gt` /
  `...expression_distribution.yaml`, matching the new `meta.id` naming
  convention above (this was the fix for a "file not found" error seen
  mid-session, twice — once when `source` was added to naming, again when
  `threshold` was added; the pattern is: any field folded into `meta.id` in
  the workflow must also be folded into the Python script's output
  filenames, since the module declares its output path as `${meta.id}.gt`).
- **`tests/test_samplesheet.csv`** — updated to the new `seeds,network,
  tissue,source,threshold` header with real values (GTEx/PAXDB/TCGA rows).

Design decisions locked in during this work (ask before changing):
- One `source`/`threshold` value applies to *all* tissues listed in a row's
  semicolon-separated `tissue` field (no per-tissue-within-a-row overrides).
- Missing `source`/`threshold` when `tissue` is set is a hard pipeline error,
  not a fallback to the global param — this was an explicit user correction
  during the session (initially planned as fallback-to-global, user rejected
  that and required the error instead).
- `--custom_filtering_file` remains a separate, global-only param — not
  made per-row in this pass.
- Doc updates (`docs/usage.md`, `nextflow_schema.json` help text) for the
  `source` column were written once, then explicitly reverted at the user's
  request ("remove the document update") — don't re-add without asking.

## 2026-08-24 → 2026-09-01 — CRAPome filtering + subworkflow extraction + Nextflow-side downloads (`tissue_specific_filtering` branch)

Since the 2026-08-10 entry above, `bin/tissue_specific_filtering.py` and
`modules/local/tissue_specific_filtering/` were renamed to
`bin/context_specific_filtering.py` / `modules/local/context_specific_filtering/`
(commit `361960b`, 2026-08-24, already pushed) — "tissue" was generalized to
"context" throughout naming since filtering sources now include non-tissue
contexts like TCGA cancer types and CRAPome. **Any reference to
`tissue_specific_filtering.py`/`.nf` elsewhere in this memory file (§1 above)
is now stale — the current name is `context_specific_filtering`.**

Three commits on `tissue_specific_filtering` as of 2026-09-01 are **unpushed**
(local-only, not yet on `origin/tissue_specific_filtering`):

- **`3214511` "implemented crapome filtering"** (2026-08-31) — adds CRAPome
  (AP-MS contaminant database, see `CLAUDE.md`) as a filtering source, run as
  a *second*, separate filtering pass rather than a `filtering_source` value
  a user picks per samplesheet row. New `CRAPOME_FILTERING` process instance
  (same `CONTEXT_SPECIFIC_FILTERING` module, aliased) runs after context
  filtering when `--filter_crapomes` is true, at a fixed
  `--crapome_filtering_threshold` (spectral-count cutoff, default 1000) —
  independent of the per-row `source`/`threshold` samplesheet columns.
  `filter_network()` in `bin/context_specific_filtering.py` special-cases
  `source == "CRAPome"`: it skips `save_expression_distribution()` (no
  meaningful expression curve for a contaminant flag) and writes a
  differently-shaped stats row (`write_crapome_filtering_statistics`) via a
  new MultiQC section `filtering_statistics_crapomes` in
  `assets/multiqc_config.yml`. The `expression_distribution` module output
  was made `optional: true` to allow this. CRAPome-filtered networks are
  `.mix()`-ed back into the main network/seeds channels alongside (or
  instead of, on top of) context-filtered ones, so a network can go through
  context filtering, CRAPome filtering, both, or neither depending on
  `params.filter_crapomes`.
- **`4ad7465` "migrated context-specific filtering into a new subworkflow"**
  (2026-08-31) — extracted the context+CRAPome filtering logic that had
  accumulated inline in `workflows/diseasemodulediscovery.nf` into a new
  `subworkflows/local/network_filtering/main.nf` (`NETWORK_FILTERING`
  workflow), keeping the top-level workflow file leaner. Same commit
  introduced a `context_key` (`context + "." + source + "." + threshold`)
  threaded through the samplesheet-parsing and CLI-input paths in
  `subworkflows/local/utils_nfcore_diseasemodulediscovery_pipeline/main.nf`
  so repeated context/source/threshold combos dedupe correctly. Also added a
  first version of `downloadFilteringFile(context, source)` — a Groovy
  helper using Nextflow's `file(url)` to fetch the GTEx/PaxDB/TCGA/CRAPome
  source files, called from within the new subworkflow.
  **Known leftover from this commit:** in `workflows/diseasemodulediscovery.nf`
  the old inline filtering block (now replaced by the `NETWORK_FILTERING(...)`
  call) was left in place but commented out with `/* ... */` rather than
  deleted — dead code that should be cleaned up.
- **`f62002e` "refactored file download to be handled by nextflow"**
  (2026-09-01) — moved *all* remaining file downloads (GTEx, PaxDB, TCGA,
  CRAPome) out of `bin/context_specific_filtering.py` (which previously used
  `requests.get` directly, including one `verify=False` call added for
  CRAPome in the prior commit to work around an SSL issue — that workaround
  is now moot since the script no longer makes HTTP requests) and into
  Nextflow itself via `downloadFilteringFile()` in the new subworkflow, using
  a `path(filtering_file)` process input. Nextflow handles retries and lets
  work-dir caching apply to downloads. `resolve_expression_file()` in the
  Python script now takes `--filtering_file` instead of a hardcoded URL
  constant per source. `downloadFilteringFile()` also gained input validation
  (`error(...)` on an unrecognized `source`, or a TCGA `context` not prefixed
  `TCGA_`) and fixed a source-name-casing mismatch (`"GTEX"/"CRAPOME"` vs. the
  actual `"GTEx"/"CRAPome"` values used elsewhere) that would have made GTEx
  and CRAPome downloads silently return `null` from the first version of the
  helper.
  **Also touches `conf/test.config`:** `seeds`/`network` test-data params
  were commented out and `validate_online`, `run_seed_perturbation`,
  `run_network_perturbation` flipped to `false`, plus a new `skip_digest =
  true`. No replacement input (e.g. a test samplesheet) was wired into
  `test.config` in this commit — **`-profile test` likely has no seeds/network
  input right now and may not run end-to-end; verify before relying on the
  test profile.**

## Open items / things to check next session

- **Push blocked:** `tissue_specific_filtering` is 3 commits ahead of
  `origin/tissue_specific_filtering` (`3214511`, `4ad7465`, `f62002e` — see
  above). Not yet pushed as of 2026-09-01.
- `conf/test.config` has `seeds`/`network` commented out with nothing
  replacing them (see `f62002e` above) — check whether `-profile test` still
  runs, or if a samplesheet-based input needs to be added there.
- Dead commented-out filtering block left in
  `workflows/diseasemodulediscovery.nf` after the `4ad7465` subworkflow
  extraction — should be deleted.
- `af736f2 "failing yaml version of network degree distribution plot"` —
  commit message implies this was left in a broken/failing state; worth
  checking if it was fixed later or still needs attention.
- A stray `.main.nf.swp` (vim swap file) was committed in both branches
  under the tissue-specific-filtering module directory — should probably be
  removed and added to `.gitignore`.
- Confirm whether `tissue_specific_filtering-multiple_filtering_sources` can
  now be deleted/archived since it's fully merged into `tissue_specific_filtering`.
- Untracked `.nf-test.log` / `.nf-test/` present in the working tree as of
  2026-09-01 (likely from a local nf-test run) — not committed, check if
  they should be `.gitignore`d.

## 2. Exploratory analysis in this repo

**2026-08-07 — split into three notebooks.** `playground.ipynb` (scratch/API
trials) and `filtering_sanity_check.ipynb` (GTEx filtering validation:
contamination analysis, expression CDF, KEGG pathway-survival heatmaps,
PaxDB BRAIN variant) were merged and re-split by purpose into:

- **`pipeline_experimentation.ipynb`** — scratch/API trials: expression
  node-filtering trials, ID-space mapping (Ensembl↔HGNC↔UniProt↔Entrez via
  pybiomart/g:Profiler/UniProt idmapping/mygene.info), ProteomicsDB API
  exploration, expression-colored graph viz trial, UCSC Xena/TCGA API trial.
- **`filtering_source_exploration.ipynb`** — analysis of the expression
  sources themselves, independent of any filtering run: GTEx expression CDF
  vs. network nodes, GTEx-vs-PaxDB tissue clustering (Spearman + cosine).
  Currently GTEx-only for CDF; PCA of expression sources is a planned
  addition, not yet implemented.
- **`filtering_results.ipynb`** — results after running the filter: GTEx +
  PaxDB filtering runs across thresholds, filtering statistics, tissue
  contamination analysis, KEGG pathway-survival heatmaps (now with an
  **"unfiltered" baseline column** showing what fraction of each pathway's
  full KEGG gene set is present in the raw network before any threshold is
  applied — added 2026-08-07 to distinguish "gene never made it into the
  network" from "gene was filtered out by the threshold"), and downstream
  disease-module Jaccard-similarity clustering.

Both original notebooks were deleted (`git rm -f`) after their content was
fully migrated; recoverable via git history if needed. `playground.ipynb`
had no prior git history before this session's commit `6be61cf`, so its
content is summarized from state rather than diffed below.

By the time of the 2026-08-07 split, `playground.ipynb` had grown well past
the 6-cell version previously summarized here (this section had gone stale —
it undercounted the file, which actually held 34 cells including a full
PaxDb/GTEx tissue-clustering analysis, ID-space mapping trials, ProteomicsDB
and UCSC Xena/TCGA API exploration, and a toy graph-coloring demo). See the
notebook bullets above for where each piece now lives; don't trust a
cell-by-cell breakdown of the old file going forward since it no longer
exists — check the three current notebooks directly instead.

`merge_paxdb.py` — standalone script (not yet part of the pipeline) that
downloads all per-tissue PaxDb integrated protein-abundance datasets for a
given organism (default `9606`, human) from `pax-db.org` (v6.1) and merges
them into one wide TSV (`id` + one column per tissue), mirroring the layout
of the GTEx `gene_median_tpm` file used by
`bin/tissue_specific_filtering.py` in the pipeline. This is the script that
produces `paxdb_all_tissues.tsv` consumed by
`filtering_source_exploration.ipynb`.

### Open items
- Decide whether `pipeline_experimentation.ipynb`'s trial-and-error cells
  (ID mapping, ProteomicsDB/Xena API exploration) are still needed now that
  the pipeline's own ID-mapping and TCGA support exist, or can be trimmed.
- PCA of expression sources (mentioned as a target analysis for
  `filtering_source_exploration.ipynb`) is not yet implemented — only the
  GTEx expression CDF and GTEx-vs-PaxDB clustering exist so far.
- `merge_paxdb.py` is currently a standalone helper in this repo; consider
  whether it should move into the pipeline fork's `bin/` alongside
  `tissue_specific_filtering.py` if PaxDb becomes a supported filtering
  source there.
- None of this session's notebook changes (§2, 2026-08-07) are committed
  yet — `playground.ipynb`/`filtering_sanity_check.ipynb` deletions are
  staged, the 3 new notebooks are untracked. Repo is also still mid-merge
  (diverged from `origin/main`, unresolved `git merge` in progress) —
  unrelated to this session's work, flagged but not touched.
