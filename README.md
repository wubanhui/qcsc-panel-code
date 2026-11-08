# Supplementary code

Reproducibility package for the manuscript

**"A qCSC-Derived Transcriptional Panel Reflects Epithelial Differentiation Rather Than Stemness in Colorectal Cancer"**

---

## 1. What is in this package

| Folder | Content |
|---|---|
| `00_config.R` | Single place where every path is declared. All scripts source this file. |
| `gene_signature/` | **Part 1** - construction of the gene signature: TCGA sample selection, differential expression, ESTIMATE adjustment, symbol and biotype filtering, the 301 -> 885 funnel summary, and the Hallmark GSEA of Fig. 1D. Scripts are numbered in the order they must be run. |
| `singlecell_clinical/` | **Part 2** - single-cell analysis, panel selection, external validation and the clinical (MMR, CMS, prognosis) associations. |
| `drug/` | **Part 3** - drug perturbation of the panel (CX-5461, Wortmannin, Emetine, Calcitriol) and the pathway-level comparison. |
| `figures/` | One script per main figure; each reads the tables written by Parts 1-3. |
| `data/` | Input data (not included; see section 4). |
| `objects/` | Intermediate `.rds` objects written by the scripts. |
| `output/` | Run logs and status tables. |
| `run_all.R` | Runs every step in order. |

---

## 2. Requirements

### Verified environment

The package was run end to end on the following combination. Use it, or a newer
one, when reproducing the results.

| Component | Version |
|---|---|
| R | **4.6.1** |
| Bioconductor | **3.23** |
| S4Vectors | **>= 0.47.6** (required by S4Arrays; 0.50.3 was tested) |
| S4Arrays | 1.12.1 |
| Seurat | 5.6.0 |
| scDblFinder | 1.26.9 |
| limma | 3.68.5 |
| clusterProfiler | 4.20.3 |
| org.Hs.eg.db | 3.23.1 |
| msigdbr | 26.1.1 |

| Package | Needed for |
|---|---|
| Seurat, SeuratObject | single-cell analysis (Part 2) |
| scDblFinder, SingleCellExperiment | doublet identification (01a, 02g) |
| limma (3.68.x) | differential expression (01d, 01f) |
| GSVA (2.6.x) | ssGSEA scoring (01c) |
| estimate | stromal and immune scores (01f) |
| clusterProfiler, org.Hs.eg.db, msigdbr | enrichment (01e, 01f, 05c, 06a, 06b) |
| survival, pROC | survival and ROC analysis (04a-04c, 05a, 05b) |
| hgu133plus2.db, hgu133a.db | Affymetrix probe annotation (04b, 04c, 07e) |
| Matrix, data.table, ggplot2, patchwork | data handling and plotting |
| ggrepel (optional) | gene labels on the volcano plot (08a); skipped when absent |

The exact versions used for the results reported in the manuscript are listed
in `sessionInfo_R4.6.1.txt`.

---

## 3. Running the pipeline

```r
# from the package root
Rscript run_all.R --dry-run          # list the steps and whether their inputs exist
Rscript run_all.R                    # run everything
Rscript run_all.R --from 03a         # start at a given step
Rscript run_all.R --only 05b         # run a single step
```

Each step runs in its own R process, writes `logs/step_<id>.log`, and the
overall status is written to `logs/run_all_status.csv`. Steps whose inputs are
absent are reported as `SKIP` rather than `FAIL`.

---

## 4. Input data

The package contains no data. All inputs are read from `DATA_ROOT`, which
defaults to `<package>/data` and can be overridden with `QCSC_DATA_ROOT`.
Individual locations can be overridden as well; the full list is at the top of
`00_config.R`.

Expected layout of `DATA_ROOT`:

```
data/
├── scRNA/            GSE132465 and GSE188711 raw matrices; QC_stage4_* objects
├── bulk/             GSE144735, GSE39582 and GSE14333 series matrices
├── lincs/            LINCS_UP_protein_coding.txt, LINCS_DOWN_protein_coding.txt
├── deg_gsva/         TCGA differential expression tables and signature lists
├── drug/             GSE201899, GSE18005, GSE5486, GSE100785 series matrices
├── MMR_external_validation/
└── gencode.v50.basic.annotation.gtf.gz
```

`intermediate/two_cohort_panel/panel_data.rds` (about 127 MB) holds the paired
two-cohort expression object. It is not shipped with this package; it is
rebuilt from the raw data by steps 02a and 02b.

`singlecell_clinical/01c_download_external_scrna.R` is a convenience script that
downloads the independent cohort GSE144735. It needs network access and is
therefore not one of the 44 steps. If you already have the series matrix, place
it under `bulk/` and skip this script.

### File formats that are accepted

The scripts adapt to the file as distributed, so either form works:

* the TCGA expression matrix may be a plain `tcga_RSEM_Hugo_norm_count` TSV or
  the gzip-compressed `tcga_RSEM_Hugo_norm_count.gz`;
* series matrices may be plain or gzip-compressed;
* gene lists are read line by line and may carry a UTF-8 byte-order mark.

### Where the results are written

Inputs are only read. Every output goes to `tables/`, `figures/`, `objects/` or
`output/` inside the package. When an input directory is writable the script
keeps the intermediate files next to the inputs, as the original analysis did;
otherwise it falls back to the package `output/` folder, so a read-only data
mount does not stop the run.

---

## 5. Order of execution

**Part 0 - Single-cell QC objects**

| Step | Script | Output |
|---|---|---|
| 00a | `01a_singlecell_qc_objects.R` | QC objects rebuilt from the raw GEO files |
| 00b | `01b_build_qc_objects.R` | per-cohort QC statistics |

**Part 1 - Gene signature**

The order below follows the original analysis: the protein-coding signature is
built first, then TCGA is matched and scored, the differential expression and
its GSEA are computed, the ESTIMATE adjustment is applied, and only then are the
gene symbols cleaned and the funnel summarised.

| Step | Script | Original script | Output |
|---|---|---|---|
| 01a | `01a_signature_protein_coding.R` | signature_construction.R | 301 candidate genes -> 290 protein-coding signature |
| 01b | `01b_tcga_sample_matching.R` | tcga_matching.R | 380 primary tumours, `TCGA_CRC_380_samples.RData` |
| 01c | `01c_tcga_scoring.R` | gsva.R | ssGSEA scores, median split 190 / 190 |
| 01d | `01d_tcga_differential.R` | deg.R | 2,795 DEG (2,292 up / 503 down) |
| 01e | `01e_hallmark_gsea.R` | gsea.R | Hallmark NES for the unadjusted model |
| 01f | `01f_estimate_adjust.R` | robustness_estimate.R | 1,595 ESTIMATE-adjusted DEG; 1,295 concordant; GSEA comparison |
| 01g | `01g_symbol_and_biotype_filter.R` | gene_cleaning_optimized.R | 885-gene consensus set (768 up / 117 down) |
| 01h | `01h_gene_funnel_summary.R` | - | `tables/gene_funnel_summary.csv` |

**Part 2 - Single cell, panel and clinical association**

| Step | Script | Output |
|---|---|---|
| 02a-02f | epithelial extraction, CNV and expression classifiers, malignancy integration, normal epithelium | cell-level objects |
| 02g | doublet rate per sample for GSE132465 (QC descriptor; no cells removed) | `GSE132465_doublet_rate_*.csv` |
| 03a | stability selection (30 stratified 70/30 splits) | 20-gene panel |
| 03b | panel concordance, fidelity, decile AUC, TCGA scores | `final_panel_stable.csv` |
| 03c | up-arm-only sensitivity, patient-level bootstrap | sensitivity tables |
| 03d | signed direction vector | `panel_20_direction_vector.csv` |
| 03e | explicit-threshold panel-versus-full discrimination (AUC, precision, recall) | `panel_discrimination_102030.csv` |
| 03f | final panel figures | figures |
| 03g | one-page summary; every number read from an existing result file | `FINAL_RESULTS_SUMMARY.csv` |
| 04a-04c | GSE144735, GSE39582, GSE14333 | external validation tables |
| 05a-05b | MMR discrimination and hardening | AUC, threshold, subgroups, leave-one-gene-out |
| 05c | panel annotation, CMS | annotation tables |
| 06a-06b | genome-wide correlation and ORA background | 25,775-gene correlation, GO tables |

**Part 3 - Drug perturbation**

| Step | Script | Output |
|---|---|---|
| 07a | `07a_CX5461.R` | CX-5461 delta score |
| 07b | `07b_Wortmannin.R` | Wortmannin delta score |
| 07c | `07c_Emetine.R` | Emetine delta score |
| 07d | `07d_Calcitriol.R` | Calcitriol tumour-versus-normal delta |
| 07e | `07e_pathway_sets.R` | pathway-level concordance |

**Figures**

| Step | Script | Figure |
|---|---|---|
| 08a | `08a_Fig1_volcano.R` | Fig. 1A |
| 08b | `08b_Fig2B_umap_lineage.R` | Fig. 2B (UMAP of all GSE188711 cells by lineage) |
| 08b2 | `08b_Fig2E_umap_zoom.R` | Fig. 2E (zoom on the epithelial clusters 5 and 14) |
| 08c | `08c_Fig3_normal_tumour.R` | Fig. 3A, 3B |
| 08d | `08d_Fig4_prognosis_MMR.R` | Fig. 4 |
| 08e | `08e_Fig5_annotation.R` | Fig. 5 |
| 08f | `08f_Fig6_drug_ROC.R` | Fig. 4E-G |

---

## 6. Reproducibility

`VALIDATION.md` reports the result of an end-to-end run and a table-by-table
comparison with the values reported in the manuscript.

---

## 7. Notes

### Gene sets used by step 07e

The pathway analysis in `07e_pathway_sets.R` uses three MSigDB collections
(2024.1.Hs): Hallmark (`h.all`), GO Biological Process (`c5.go.bp`) and
canonical pathways (`c2.cp`). Place the three files in
`<QCSC_DATA_ROOT>/mechanism_standard_R/gmt/` as `h.all.gmt`, `c5.go.bp.gmt` and
`c2.cp.gmt`; the script looks there first and only downloads them from
https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/ when they
are absent. Steps 07a-07d do not need any gene set file.


### Run-to-run variability of the GSEA enrichment scores

Every step is seeded (`set.seed(1234)`, and the 30 stability-selection splits use
seeds 1001-1030), so all counts, correlations, AUCs, hazard ratios, P values and
gene-set mean scores are reproduced exactly.

The one exception is the normalised enrichment score (NES) produced by
`clusterProfiler::GSEA` in steps 01e and 01f. Roughly 5% of the ranked genes carry
tied test statistics, and the order of tied genes changes the permuted
enrichment score. Running the same GSEA twice on the same ranked list gives an
identical NES, but two independent pipeline runs can differ by up to about 0.03.
The direction of every score and the overall pattern are unaffected: in all runs
49 of the 50 Hallmark pathways change in the same direction between the
unadjusted and the ESTIMATE-adjusted model, and
HALLMARK_OXIDATIVE_PHOSPHORYLATION is the only pathway that flips. Its adjusted
NES is about +1.17 in the delivered code (the manuscript quotes +1.14 from the
original analysis, a difference of 0.03).

* Every script locates the package root automatically, so it can be started
  from any working directory.
* No absolute path appears anywhere in the scripts.
* Random seeds: `set.seed(1234)`; the 30 stability-selection splits use seeds
  1001-1030.
* The panel is deterministic: the same 20 genes and the same selection
  frequencies are obtained on every run.

---

## 8. Licence

MIT. See `LICENSE`.
