# I. Gene expression analyses

To quantify gene expression across cells we developed an in-house pipeline. This pipeline can be used to map reads from VASA-seq, 10x and SMART-seq data. The full description of this pipeline can be found at the [a_Mapping](a_Mapping) folder. Main difference between this and other pipelines is that it recovers multi-mapper reads falling into gene bodies using a hierarchical strategy, it annotates unspliced and spliced transcripts separatedly, and it keeps biotype identity in the gene name. 

Next, we combined the use of scrublet and scanpy with in-house code to perform downstream single-cell RNA sequencing analysis (cell filtering, doublet detection, cell cycle regression, umap projection, cell type calling, differential gene expression analysis, etc). These scripts can be found at [b_Analysis](b_Analysis).
