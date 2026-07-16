# Overview

In addition to the transcriptome assembly discussed, we ran [MicroExonator](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02246-2) to identify short unannotated exons. Install instructions and general information can be found at MicroExonator's [documentation](https://microexonator.readthedocs.io/en/latest/) page. Since MicroExonator is still in active development, you can access the specific version we used for this project by accessing a branch called `joe` with the following command:

```{bash}
git checkout joe
```
Current versions may produce similar results, but the configuration files needed will be slighly different to the one provided here.

# Analysis

Similarly to the transcriptome assembly analyses, MicroExonator relies on several configuration and data files to run. We ran the [Discovery](https://microexonator.readthedocs.io/en/latest/discovery_and_quantification.html) module using the files that can be found at the [Build](/Build) folder. The following command was used to execute the analysis:

```{bash}
snakemake -s MicroExonator.skm --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" --use-conda -k -j 2000 Report/out.high_quality.txt
```



