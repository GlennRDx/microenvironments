## How to run through analysis:
1. Clone repo
2. Download janesick_2023 directory from google drive and place in 01_data/
3. Run through analysis files sequentially located in 02_analysis/


## Directory and files tree structure:
```bash
microenvironments. 
├── 01_data
│   └── janesick_2023
│       ├── processed
│       │   ├── fit_lda.RDS
│       │   ├── genes_to_keep.RDS
│       │   ├── iss_meta.RDS
│       │   ├── iss_obj_invasive.RDS
│       │   ├── shared_genes.RDS
│       │   ├── train_genes.RDS
│       │   ├── wts_obj.RDS
│       │   └── wts_obj_invasive.RDS
│       ├── raw
│       │   ├── Fig2a_scFFPE-seq_UMAP.csv
│       │   ├── Fig3e-j_Xenium.csv
│       │   ├── Fig3k_Heatmap.csv
│       │   ├── ISS_cell_types.txt
│       │   ├── ISS_cells.csv
│       │   ├── ISS_gene_expression.txt
│       │   ├── scRNA_cell_types.txt
│       │   ├── scRNA_gene_expression.txt
│       │   └── scRNAseq_data_quantile_matched.txt
│       └── raw_clean
│           ├── iss_expr_clean.RDS
│           ├── iss_meta_clean.RDS
│           ├── wts_expr_clean.RDS
│           ├── wts_meta_clean.RDS
│           └── wts_obj_clean.RDS
├── 02_analysis
│   ├── 0.0.1_packages.Rmd
│   ├── 0.0.2_data_cleaning.Rmd
│   ├── 0.1.1_scRNA_seq analysis.Rmd
│   ├── 0.2.1_ISS surface_core_annotation.Rmd
│   ├── 0.3.1_lda_surface_core_classification.Rmd
│   ├── 0.4.1_whole_transcriptome_surface_prediction.Rmd
│   ├── 0.5.1_TME_prediction_validation.Rmd
│   ├── 0.6.1_GO enrichment.Rmd
│   ├── Archive
│   │   ├── 0.3.2_svr_model.Rmd
│   │   ├── 0.3.3_random_forest_test.Rmd
│   │   └── 0.5.2_TME_prediction_validation_DESeq2.Rmd
│   └── config.R
├── 03_figures
├── README.md
└── microenvironments.Rproj
```