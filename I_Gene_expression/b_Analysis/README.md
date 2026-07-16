# VASAseq Analysis Pipeline

The same pipeline was used to analyze the data from (Pijuan-Sala et al, Nature 2018)[https://www.nature.com/articles/s41586-019-0933-9], changing only the relevant parameters in the filtering steps. 

## Contents

### Scripts
- `01_qc_checks_HPCversion.py`: In this notebook we perform the quality checks (QC) of the VASA libraries (mouse embryo data: E6.5, E7.5, E8.5 and E9.5) and identify doublet cells with `scrublet`. The notebook takes as an input the merged tables (in feather format) for both counts obtained using reads mapping in the whole gene bodies, or counts obtained from reads that fall into the 80% side of the 3' end of the gene. 

- `02_scanpy_QCxBiotype.py`: Creates 3 different scanpy objects: 1 for only spliced counts, another for only unspliced counts, and the 3rd for both spliced and unspliced counts. In the 3rd object, spliced and unspliced counts from the same gene are given different gene identifiers. Each scanpy object is further split into new objects, in which only counts from different biotypes are used (protein-coding, smallRNA, lncRNA, tRNA). For all objects, only genes present in at least 2 cells are kept. Extra QC plots showing histone content or fraction of biotypes per cell are produced. 

- `03_cellCycle.py`: Here, we identify cells in S-phase using histone content for all timepoints, separatedly and combined, based on histone expression (using only spliced read counts). We use only count tables produced with all gene counts (falling into any position along the gene body).  We further filter cells based on the parameters in `filterParams.py`. The script performs differential gene expression analysis between S-phase and non-S-phase cells, for both spliced and unspliced counts (separately) using the `scanpy` function.

- `04_scanpy_Filtered_suUMAP_Hist.py`: Performs basic scRNA-seq analysis; filtering, normalization, UMAP, differential gene expression, using mainly the steps described in scanpy.  It does this for spliced, unspliced and unspliced+spliced data simultaneously, and for each of the different biotypes ('All','ProteinCoding','lncRNA','smallRNA','TF'). An extra filtering step is included in while building the kNN graph to filter out cells that are very far from their first neighbor. As a consequence, not all cells have the same number of neighbors. This step allows us to remove extra doublets and artifacts that were missed by the previous QC steps. Filtering parameters etc can be found in the `filterParams` self-made package. User-made functions are found in the self-made `plot_aautils` and `sc_aautils` packages.

- `05_PJvsVASA.py`: First step to compare the 10x-based mouse embryo development atlas with the VASAseq-based mouse embryo atlas. The script takes as an input all the 10X and VASA-seq h5ad.gz objects produced in step 04, and finds nearest neighbors between the two atlanses in order to identify groups of cells that are equivalent. To do so, it works with both the h5ad object obtained using all reads falling in any part of the gene body, or only reads that fall in the 3'UTR. 

- `06_masterUMAP.py`: Performs projection of the full dataset (all the different timepoints) into the same UMAP. For this, cells from E6.5 are connected with cells from E6.5, while cells from E7.5 can be connected with cells from E6.5 and cells from E7.5, and cells from E8.5 can be connected only with cells with E7.5 and E8.5, etc. To connect cells between sequential timepoints, all cells are projected to the PCA space obtained with the initial timepoint. 

- `06_masterUMAP_noHist.py`: Same as before, after regressing out S-phase from the transcriptome data. 

- `07_RNAvelocity.py`: Performs RNAvelocity analysis using `scvelo` for each timepoint. 

- `08_RNAvelocity_merged.py`: Performs RNAvelocity analysis using `scvelo` on the whole dataset, and zooms in into the blood and the endothelium lineages. 

### In house python libraries
- `filterParams.py`: For each timepoint, this contains the parameters used to filter out cells with not enough sequencing quality, the number of PCA used to build the kNN graph, the resolution set to perform leiden clustering, and the dispersion selected to perform UMAP projection. These parameters are common for all libraries in each timepoint. In E8.5, one full library is excluded. As an output, several QC plots are generated. 
- `plot_aautils.py`: Collection of functions used a lot to generate graphs with our cherry-picked properties. 
- `sc_aautils.py`: Collection of functions used recurrently to analyze single cell data and that are not included in the `scanpy` package or that are wrappers of functions that are indeed included in `scanpy`.

### Input parameters
- `HistoneGenes.tsv`: list of histone genes
- `Mus_musculus_TF.txt`: List of annotated murine transcription factors.

## Detailed description of the pipeline

All the scripts can be run direclty in a LINUX/UNIX terminal using a virtual environment with python3 and all the required libraries (`pandas, numpy, collections, scanpy, scrublet, ...`) or submitting the jobs via SLURM or SGE. In some cases, a lot of memory is needed and the script takes some time, therefore job submission is recomended. 

### How to run the scripts in the LINUX/UNIX terminal

- `timepoint` indicates which timepoint are we analyzing and as an input of our scripts we set it as: E65, E75, E85, E95
- `genecoverage` indicates whether we use all reads or only reads mapping the the 3' UTR of detected genes. Can have the values `all` or `high`. We recommend to run 02 with both genecoverage options before step 03, since both outputs are required for step 05. 

```
timepoint=E65
genecoverage=all
01_qc_checks_HPCversion.py $timepoint 
02_scanpy_QCxBiotype.py $timepoint all
02_scanpy_QCxBiotype.py $timepoint high
03_cellCycle.py
04_scanpy_Filtered_suUMAP_Hist.py $timepoint $genecoverage
05_PJvsVASA.py $timepoint
07_RNAvelocity.py $timepoint
```
Once the previous steps have run for all the timepoints, we can proceed with the next scripts: 
```
06_masterUMAP.py
06_masterUMAP_noHist.py
08_RNAvelocity_merged.py
```

