# CLAUDE.md — PPI Network Filtering for Disease Module Detection

**Thesis project at TUM School of Life Sciences, Chair of Data Science in Systems Biomedicine**
**Context:** REPO4EU project (https://repo4.eu/)

---

## Project Goal

Integrate tissue- and disease-specific filtering of protein-protein interaction (PPI) networks into the existing `nf-core/diseasemodulediscovery` Nextflow pipeline to improve the biological relevance and reliability of discovered disease modules. Analyze the effect of filtering on module quality using the pipeline's built-in evaluation framework.

---

## Repository

**Upstream (nf-core):** https://github.com/nf-core/diseasemodulediscovery
**Working fork:** https://github.com/Schansiate/diseasemodulediscovery — this is the remote origin of the local `diseasemodulediscovery` repository; all practical implementation work lives here
**Feature branches:** branches prefixed with `tissue_specific_filtering` contain the GTEx filtering implementation
**Branch to work on:** `dev`
**Pipeline docs:** https://nf-co.re/diseasemodulediscovery/dev/
**Preprint (cite this):** Kersting et al., bioRxiv 2025. doi:10.1101/2025.11.20.687681

---

## Theoretical Background

### Disease Module Identification

Most human diseases arise from disruptions of specific molecular subnetworks rather than individual genes. Disease module identification (also called active module identification, AMI) formalizes this by:

1. Starting from a set of known disease-associated genes/proteins (**seeds**), sourced from databases like DisGeNET, ClinVar, or GWAS results.
2. Expanding these seeds through a background PPI network (the **interactome**) into a connected subnetwork — the **disease module** — that reflects the molecular context of the disease.
3. Evaluating the biological plausibility of the resulting module.

Algorithms implemented in the pipeline:
- **DOMINO** — active module identification with reduced false-call rate (Levi et al., Mol Syst Biol 2021)
- **DIAMOnD** — iterative connectivity-based seed expansion (Ghiassian et al., PLoS Comput Biol 2015)
- **ROBUST / ROBUST bias-aware** — prize-collecting Steiner trees (Bernett et al., Bioinformatics 2022)
- **1st Neighbors** — direct neighborhood of seed genes
- **Random Walk with Restart (RWR)**

### Why PPI Network Quality Matters

PPI networks are the backbone of module detection, but they have well-documented issues:

**Problem 1 — Tissue/disease agnosticism:** A standard PPI network lists all known interactions regardless of whether both proteins are expressed in the tissue or cell type relevant to the disease. If a protein is not expressed in liver tissue, it cannot participate in liver-specific interactions — yet it appears in the network and can distort module shape.

**Problem 2 — Study bias:** PPIs are not uniformly characterized; highly-studied "hub" proteins accumulate more interaction evidence simply because they are studied more. This creates artificial hubs and a scale-free degree distribution that biases module detection (Blumenthal et al., eLife 2024).

**Problem 3 — Experimental false positives:** Affinity purification–mass spectrometry (AP-MS) experiments that identify PPIs are contaminated by non-specific background binders. CRAPome catalogs ~300 common contaminants from these experiments (Mellacheruvu et al., Nat Methods 2013).

### Tissue-Specific Filtering — Theory

The core principle: a PPI can only be biologically active if both interacting partners are expressed in the tissue or condition of interest. Therefore, removing unexpressed proteins (and all their edges) from a PPI network creates a **context-specific subnetwork** that is more appropriate for studying tissue-relevant disease mechanisms.

Evidence for this approach:
- PPIs are themselves tissue-specific (Laman Trip et al., Nat Biotechnol 2025)
- Disease modules are tissue-specific (Kitsak et al., Sci Rep 2016)
- IID demonstrates improved disease-relevance using tissue-filtered networks (Kotlyar et al., Nucleic Acids Res 2016)

**Filtering pipeline (conceptual):**
```
Raw PPI network
       │
       ▼
Expression database (GTEx / Recount3 / PaxDB / ProteomicsDB)
       │  select tissue/condition
       ▼
List of expressed genes/proteins (above threshold)
       │
       ▼
Remove nodes NOT in expressed set  ──→  Remove all incident edges
       │
       ▼
Context-specific PPI network  ──→  Feed into disease module pipeline
```

**Threshold choice:** Typically a gene is considered expressed if its median TPM (or RPKM) across samples of the tissue exceeds a cutoff (commonly TPM > 0.5 or TPM > 1). For protein abundance (PaxDB/ProteomicsDB), a ppm abundance threshold can be used. The exact threshold is a parameter.

---

## The Nextflow Pipeline

### What it does

`nf-core/diseasemodulediscovery` is an all-in-one Nextflow pipeline that:
- Takes seed gene/protein files and a PPI network as input
- Runs multiple AMI tools in parallel
- Evaluates modules for topology, functional coherence, seed rediscovery, and robustness
- Annotates results with drug/disease information from NeDRexDB
- Exports modules to BioPAX format and the Drugst.One visualization tool
- Generates a summary MultiQC report

### Quick Start

```bash
# Install Nextflow (Java 11+ required)
curl -s https://get.nextflow.io | bash

# Test the pipeline
nextflow run nf-core/diseasemodulediscovery \
  -r dev \
  -profile docker,test \
  --outdir ./test_results

# Run with your own data
nextflow run nf-core/diseasemodulediscovery \
  -r dev \
  -profile docker \
  --seeds ./seeds.txt \
  --network string_min900 \
  --id_space entrez \
  --outdir ./results
```

> **Note:** Always use `-r dev` as the pipeline is still under active development and has not had a stable release.

### Input Formats

**Seeds file** (`--seeds`): plain text, one gene/protein ID per line, no header.
```
2717
175
4669
```

**Network file** (`--network`): CSV edge list (no header), or `.gt` / `.graphml` / `.dot`, or a built-in keyword (see table below).
```
3920,5476
113457,4214
```

**ID spaces** (`--id_space`): `entrez`, `ensembl`, `hgnc`, or `uniprot`.

### Built-in PPI Networks

| Keyword | Source | Nodes | Edges | Notes |
|---|---|---|---|---|
| `string_min900` | STRING v12.0 | 11,971 | 93,559 | Score > 0.9 |
| `string_min700` | STRING v12.0 | 15,788 | 224,045 | Score > 0.7 |
| `string_physical_min900` | STRING v12.0 | 7,722 | 34,141 | Physical only, score > 0.9 |
| `string_physical_min700` | STRING v12.0 | 10,465 | 78,878 | Physical only, score > 0.7 |
| `biogrid` | BioGRID 4.4.242 | 18,101 | 865,553 | |
| `hippie_high_confidence` | HIPPIE v2.3 | 13,246 | 112,202 | Score > 0.73 |
| `hippie_medium_confidence` | HIPPIE v2.3 | 16,613 | 637,499 | Score > 0.63 |
| `iid` | IID 18.03.2025 | 19,598 | 1,202,716 | |
| `nedrex` | NeDRexDB 18.03.2025 | 18,718 | 935,139 | Experimentally validated |
| `nedrex_high_confidence` | NeDRexDB 18.03.2025 | 12,827 | 95,944 | Score > 13.5 |

Networks are hosted on Zenodo (https://zenodo.org/records/15049754) and downloaded automatically.

### Key Parameters

| Parameter | Description |
|---|---|
| `--seeds` | Path to seed file(s), comma-separated for multiple |
| `--network` | Network file(s) or keyword(s), comma-separated |
| `--id_space` | ID space: `entrez`, `ensembl`, `hgnc`, `uniprot` |
| `--outdir` | Output directory |
| `--run_seed_perturbation` | Enable seed perturbation robustness analysis |
| `--run_network_perturbation` | Enable network perturbation robustness analysis |
| `--skip_domino` | Skip DOMINO |
| `--skip_diamond` | Skip DIAMOnD |
| `--skip_robust` | Skip ROBUST |

Use `-params-file params.yaml` to avoid long CLI commands. Do NOT use `-c` to pass parameters.

### Pipeline Architecture (nf-core)

The pipeline follows nf-core conventions:
```
main.nf                  # Entry point
workflows/               # Top-level workflow logic
subworkflows/            # Reusable sub-workflows
modules/                 # Individual process definitions
conf/                    # Configuration files
assets/                  # MultiQC config, schema, etc.
bin/                     # Helper scripts
```

Each process in `modules/` specifies its own container. Containers are pulled automatically with `-profile docker` or `-profile singularity`.

### Nextflow Concepts Relevant to This Project

**Process:** A single computation unit. Defined with `process { ... }`. Runs in its own work directory.

**Channel:** Data flowing between processes. Use `Channel.fromPath()` for files, `Channel.value()` for scalars.

**Workflow:** Connects processes with channels. The main logic lives in `workflows/diseasemodulediscovery.nf`.

**Profile:** A named collection of config settings. `-profile docker` uses Docker containers; `-profile singularity` uses Singularity.

**Resume:** `-resume` skips already-completed processes using cached results. Essential for iterative development.

**Params file:** A YAML/JSON file with pipeline parameters passed via `-params-file`. Preferred over long CLI flags.

**MultiQC integration:** Results from individual processes are aggregated by MultiQC into an HTML report. Integration requires writing output to the MultiQC input channel.

---

## Expression Databases for PPI Filtering

### GTEx (primary target for this project)

**What it is:** The Genotype-Tissue Expression project. Bulk RNA-seq data from ~54 human tissues collected from post-mortem donors. The canonical reference for human tissue-specific gene expression.

**URL:** https://gtexportal.org/

**Key data:**
- Gene-level TPM expression matrices per tissue
- ~17,000 samples across 54 tissues (v10)
- Provides median TPM per gene per tissue — ideal for filtering

**Access:**
- Portal download (requires registration for raw data, but processed matrices are open)
- GTEx API: https://gtexportal.org/api/v2 — can query by tissue and gene
- Bulk download via Google Cloud Storage or AWS

**Relevant IDs:** Gene expressions are annotated with Ensembl gene IDs (ENSG*). Will require mapping to the ID space used by the pipeline (Entrez, HGNC, UniProt).

**Integration plan:**
- Pipeline parameter `--gtex_tissue` (e.g., `"Liver"`, `"Brain - Frontal Cortex (BA9)"`)
- Download or cache the median-TPM matrix for the selected tissue
- Filter out all genes with median TPM below a configurable threshold (`--expression_threshold`, default 0.5 or 1)
- Map remaining expressed genes to the pipeline's ID space
- Remove unexpressed nodes (and their edges) from the PPI before module detection

**MultiQC reporting:** Report number of nodes before/after filtering, percentage of seeds retained.

### Recount3

**What it is:** A large harmonized collection of human and mouse RNA-seq studies processed through a uniform pipeline. Covers GTEx, TCGA, and thousands of SRA studies.

**URL:** https://rna.recount.bio/ / R package: `recount3`

**Advantage over GTEx:** Broader tissue/disease/condition coverage (including disease states); useful when GTEx doesn't have the specific condition needed.

**ID space:** Ensembl gene IDs.

### PaxDB (Protein Abundance Database)

**URL:** https://pax-db.org/

**What it is:** Curated protein abundance data (in parts-per-million, ppm) from proteomics experiments across organisms and tissues. v6.0 (2025) uses LLM-assisted curation.

**Advantage:** Uses actual protein abundance rather than mRNA, which better captures translational regulation and post-transcriptional effects.

**ID space:** STRING IDs / UniProt. Requires mapping.

### ProteomicsDB

**URL:** https://www.proteomicsdb.org/

**What it is:** Human proteome data integrated from many MS-based proteomics experiments. Includes cell line and tissue data. REST API available.

**Advantage:** Broad tissue coverage with protein-level quantification.

**ID space:** UniProt ACs.

### CRAPome (False Positive Filtering)

**URL:** https://reprint-apms.org/ (replaces original crapome.org)

**What it is:** A repository of AP-MS contaminant proteins. Contains ~300 proteins frequently detected as non-specific background in AP-MS experiments — these interactions are likely false positives.

**Usage:** Remove edges where one or both partners are high-frequency contaminants (appears in > N% of control AP-MS experiments).

### TCGA (via UCSC Xena)

**URL:** https://xenabrowser.net/datapages/

**What it is:** The Cancer Genome Atlas — expression data from tumor samples. Useful for filtering PPI networks in cancer-specific disease module analyses.

**Access:** UCSC Xena Browser provides bulk download of TCGA RNA-seq matrices.

---

## Evaluation Framework

The pipeline's built-in evaluation covers:

| Metric | Tool/Method | What it tests |
|---|---|---|
| Over-representation analysis | g:Profiler | Are module genes enriched in known pathways? |
| Functional coherence | DIGEST | Do module genes share function? |
| Network topology | graph-tool | Module size, density, connectivity |
| Module overlaps | custom | Consistency across AMI tools |
| Seed rediscovery | perturbation | Does the module recover held-out seeds? |
| Network robustness | perturbation | Is module stable under edge removal? |

**For this project, additionally evaluate:**
- Does tissue-specific filtering improve pathway relevance for tissue-specific diseases?
- Are tissue-specific pathways more enriched in filtered vs. unfiltered modules?
- What fraction of seed genes are retained after filtering? (seeds lost = seeds in unexpressed genes)
- How do module sizes change with filtering?

---

## Suggested Test Diseases

Choose diseases with clear tissue specificity vs. systemic/broad diseases for comparison:

**Tissue-specific (filtering should help):**
- Non-alcoholic fatty liver disease (NAFLD) → liver
- Alzheimer's disease → brain
- Type 2 diabetes → pancreas/adipose
- Breast cancer → breast tissue
- Idiopathic pulmonary fibrosis → lung

**Systemic (filtering effect less clear):**
- Rheumatoid arthritis
- Systemic lupus erythematosus
- Type 1 diabetes

Seed genes can be obtained from DisGeNET (https://www.disgenet.org/), ClinVar, or GWAS Catalog.

---

## Development Guidelines

### Filtering implementation

All filtering (GTEx and other expression databases) is handled by the `TISSUE_SPECIFIC_FILTERING` process in the working fork. New data resources should be integrated into this process rather than adding separate pipeline steps.

### ID space mapping

The pipeline handles multiple ID spaces (Entrez, Ensembl, HGNC, UniProt). GTEx uses Ensembl gene IDs. Mapping resources:
- `biomaRt` R package for Ensembl ↔ Entrez ↔ HGNC
- UniProt ID mapping API: https://rest.uniprot.org/idmapping/
- mygene.info API for gene ID mapping

Ensure the mapping step accounts for many-to-many relationships and logs unmapped IDs.

---

## Key References

| Topic | Citation |
|---|---|
| Pipeline paper | Kersting et al., bioRxiv 2025 (doi:10.1101/2025.11.20.687681) |
| Nextflow | Di Tommaso et al., Nat Biotechnol 2017 (doi:10.1038/nbt.3820) |
| nf-core | Ewels et al., Nat Biotechnol 2020 (doi:10.1038/s41587-020-0439-x) |
| DOMINO | Levi et al., Mol Syst Biol 2021 (doi:10.15252/msb.20209593) |
| DIAMOnD | Ghiassian et al., PLoS Comput Biol 2015 (doi:10.1371/journal.pcbi.1004120) |
| ROBUST | Bernett et al., Bioinformatics 2022 (doi:10.1093/bioinformatics/btab876) |
| GTEx | GTEx Consortium, Nat Genet 2013 (doi:10.1038/ng.2653) |
| Recount3 | Wilks et al., Genome Biol 2021 (doi:10.1186/s13059-021-02533-6) |
| PaxDB v6.0 | Huang et al., Nucleic Acids Res 2025 (doi:10.1093/nar/gkaf1066) |
| ProteomicsDB | Picciani et al., Nucleic Acids Res 2025 (doi:10.1093/nar/gkaf1265) |
| CRAPome | Mellacheruvu et al., Nat Methods 2013 (doi:10.1038/nmeth.2557) |
| Tissue-specific PPIs | Laman Trip et al., Nat Biotechnol 2025 (doi:10.1038/s41587-025-02659-z) |
| Tissue-specific disease modules | Kitsak et al., Sci Rep 2016 (doi:10.1038/srep35241) |
| IID (tissue-filtered PPI) | Kotlyar et al., Nucleic Acids Res 2016 (doi:10.1093/nar/gkv1115) |
| Study bias in PPIs | Blumenthal et al., eLife 2024 (doi:10.7554/eLife.99951) |
| MultiQC | Ewels et al., Bioinformatics 2016 (doi:10.1093/bioinformatics/btw354) |
| STRING | Szklarczyk et al., Nucleic Acids Res 2023 (doi:10.1093/nar/gkac1000) |
| BioGRID | Oughtred et al., Protein Sci 2021 (doi:10.1002/pro.3978) |
| HIPPIE | Alanis-Lobato et al., Nucleic Acids Res 2017 (doi:10.1093/nar/gkw985) |
| NeDRex | Sadegh et al., Nat Commun 2021 (doi:10.1038/s41467-021-27138-2) |

---

## Project Timeline

### Practical Work — COMPLETED
- [x] Read literature listed in references above
- [x] Evaluate GTEx API and download options
- [x] Understand GTEx data structure (tissues, sample IDs, TPM matrices, identifier scheme)
- [x] Implement GTEx PPI filtering module in Nextflow (see `tissue_specific_filtering` branches in the working fork)
  - [x] Software container
  - [x] `--gtex_tissue` pipeline parameter
  - [x] MultiQC stats (nodes before/after, seeds retained)
  - [x] ID space mapping
  - [x] Sample sheet integration
- [x] Poster presentation

**Initial test runs** with GTEx data (Alzheimer's disease) are available in `alzheimer_symbol_runs/` in this repo.

### Thesis Work (~4 months)
- [ ] Select test diseases (tissue-specific vs. systemic)
- [ ] Evaluate filtering impact using pipeline evaluation metrics
  - module overlap, seed rediscovery, ORA, tissue-specific pathway enrichment
- [ ] Extend pipeline
  - [x] Cumulative expression distribution in MultiQC
  - [x] Color expression values in graph visualizations
  - [x] Custom user input (column median for filtering)
- [ ] Evaluate additional databases: Recount3, PaxDB, ProteomicsDB, CRAPome, CoBiNet, TCGA
- [ ] Integrate additional filtering approaches (time-dependent)
- [ ] Extend evaluation framework
- [ ] Write thesis

---

## Useful Links

- Pipeline GitHub: https://github.com/nf-core/diseasemodulediscovery
- Pipeline docs: https://nf-co.re/diseasemodulediscovery/dev/
- Nextflow docs: https://www.nextflow.io/docs/latest/
- nf-core docs: https://nf-co.re/docs/
- nf-core Slack `#diseasemodulediscovery`: https://nf-co.re/join/slack
- GTEx portal: https://gtexportal.org/
- GTEx API: https://gtexportal.org/api/v2
- Network preparation repo: https://github.com/REPO4EU/network_preparation
- REPO4EU project: https://repo4.eu/
- CoBiNet (non-interactors): https://www.cobinet.ai/
- UCSC Xena (TCGA): https://xenabrowser.net/datapages/