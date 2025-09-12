
## Single cell analysis pipeline

## Suggested filestructure

```
/sc_rna/
├─ H4_A/
│  └─ H4_A_CogentDS_Analysis_Processed.rds
├─ H4_C/
│  └─ H4_C_CogentDS_Analysis_Processed.rds
├─ H4_D/
│  └─ H4_D_CogentDS_Analysis_Processed.rds
├─ H8_A/
│  └─ H8_A_CogentDS_analysis_Processed.rds
├─ H8_C/
│  └─ H8_C_CogentDS_analysis_Processed.rds
├─ H8_D/
│  └─ H8_D_CogentDS_analysis_Processed.rds
├─ H4_HEK/
│  └─ H4_HEK_CogentDS_Analysis_Processed.rds
├─ H8_HEK/
│  └─ H8_HEK_CogentDS_Analysis_Processed.rds
│
├─ scripts/
│  └─ scrna_ipsc/
│     ├─ scrna_pipeline.R          # your main pipeline
│     ├─ setup_pipeline.R          # package installer
│- refs/
│  └─ olah_2020_microglia.h5ad  # reference
│
├─ results/                         # auto-created by the pipeline
│  ├─ logs/                         # pipeline_YYYYMMDD_HHMMSS.log
│  ├─ deseq2/                       # pseudobulk DE tables (two variants)
│  ├─ DE_cluster/                   # per-cluster DE CSVs
│  ├─ gsea/                         # Hallmark fgsea CSVs (two variants)
│  ├─ variance/                     # variancePartition outputs
│  ├─ composition/                  # Dirichlet regression summary
│  ├─ HEK_QC/                       # HEK cross-chip correlation CSV
│  ├─ labels/                       # transferred labels + proportions
│  └─ markers/                      # glmnet minimal panel
│
└─ figures/                         # auto-created by the pipeline
   ├─ UMAP/                         # umap_*.png, comp-labeled UMAPs
   ├─ DE/                           # MA_*.png, Volcano_labeled_*.png, top genes by cluster
   ├─ GSEA/                         # hallmark_dotplot_*.png (two variants)
   ├─ GO/
   │  ├─ variant_labelMG/
   │  │  ├─ H4_top100/              # dotplot_BP.png, treemap_or_bar_BP.png, emap_BP.png, go_table_BP.csv
   │  │  ├─ H8_top100/
   │  │  ├─ DE_up_in_H8/
   │  │  └─ DE_down_in_H8/
   │  └─ variant_MGclusters/
   │     ├─ H4_top100/
   │     ├─ H8_top100/
   │     ├─ DE_up_in_H8/
   │     └─ DE_down_in_H8/
   ├─ HEK_QC/                       # HEK plots
   ├─ Celltypes/                    # composition stacks
   ├─ Composition/                  # subset change lines
   └─ markers/                      # ROC plot
```
