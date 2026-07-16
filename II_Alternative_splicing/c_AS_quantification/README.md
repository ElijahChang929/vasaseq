# Overview

To quantify alterantive splicing events we used [Whippet](https://github.com/timbitz/Whippet.jl) through the downstream module provided by MicroExonator. For additional details of MicroExonator's scRNA-seq module see its [documentation](https://microexonator.readthedocs.io/en/latest/single_cell_analysis.html) page.

# Analysis

We use MicroExonator's scRNA-seq module to perform pairwise comparisons between cell clusters and to obtain a PSI matrix for further downstream analyses.

## Pairwise comparisons

Neccesary files for this module can be found in the [Build](Build/Pairwise_comparisons) folder. We used the following command to perform these comparisons:

```{bash}
snakemake -s MicroExonator.skm --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" --use-conda -k -j 2000 snakepool
```

## PSI matrix

The same files as for the previous module are required, except the `config.yaml` which can be found at the [Build](Build/PSI_matrix) folder. We used the following command to obtain this PSI matrix:

```{bash}
snakemake -s MicroExonator.skm --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" --use-conda -k -j 2000 merge_quant_by_cluster_PSI
```

