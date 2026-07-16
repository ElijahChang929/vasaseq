# Overview

To perform transcriptome assembly using VASA-seq data we extended some of the workflows that were initally developed for demonstrative purposes in the Hemberg group. The [repository](https://github.com/hemberg-lab/RNA_seq_snakepipes/) of these workflows is available on GitHub, under a branch called `joe`.

This and all the computational workflows used in the alternative splicing analyses were implemented using Snakemake<sup>[1](https://snakemake.readthedocs.io/en/stable/)</sup>, a popular workflow manager that enables the creation of scalable and repoducible computational workflows.  


# Obtaining the code

To obtain the code from our repository use:

```{bash}
git clone https://github.com/hemberg-lab/RNA_seq_snakepipes/
```

Then to access the right branch, you should run:

```{bash}
git checkout joe
```

# Setup

To reproduce these analyses, snakemake needs to be installed (see  [installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)). Moreover, specific files needs to be plased inside the cloned folder, which are described bellow. The complete set of files used to run the analysis of VASA-seq data can be found at the [Build](Build/) folder.


## config.yaml

This YAML file contains the information that is required for configuring the pipeline. This file needs to be placed inside the clone's folder and the following file paths needs to be defined:


- Genome : Here you should include the path to genome (multi-fasta file containing all chromosomes)
- GTF : Path to a GTF file  - Please check that both Genome and GTF have the same chromosome names. For example European names are 1, 2, 3 …while American are Chr1, Chr2, Chr3 …


Other relevant variables:

- samples : samples.tsv
- units : units.tsv
- group_by_cluster : True


## units.tsv

The specific name and path of this file should match with the one defined in `cluster.yaml`. The folowing file should have tab-separated values (TSV) that should contain the following named columns:


- sample : ID of each cell / barcode
- unit : as we are processing single end reads, this should be 1. I have not developed this pipeline for paired ends, but I have included this parameter  if needed in the future. Leave all the values as 1.
- path : path to each fastq.gz file. We accept only compressed gz files, as it is the most efficient way to store the data and read it with reasonable speed. Your files will only be read by the pipeline and temporary files will be created inside the RNA_seq_snakepipes folder. It is advisable to always keep additional copies of the original data.

## samples.tsv

Similarly to `units.tsv`, this also a TSV file that should contain the following columns

- sample : ID of each cell / barcode
- pools : Number of sub-divisions per cluster. In our case, we did not subdivide the clusters as they were already corresponding to over 40 divisions of the whole data. The concept is to boost the coverage for the initial assembly steps. So merging all cells from the same cluster made sense for us, but if you have bigger datasets than the one used in the paper, maybe you could try to generate several pseudo pools per cluster by increasing this variable. This number must be the same for all cells belonging to the same cluster
- condition : In our case this should be the cluster ID. Leiden clusters worked fine for us.


## cluster.json

This is a file that snakemakes needs to communicate with your HPC queuing system. This is probably one of the steps that will take some time to figure out how to adapt for your case, unless you know people who are already using snakemake in your cluster. In our case our queuing system was LSF. Inside this file you will see the CPU and memory assignment for the rules, with some default and custom definitions for certain rules. For more information the official snakemake documentation on this topic:

https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration-deprecated



# Run 

To execute this pipeline, activate snakemake's conda enviroment and provide the appropiate cluster variables and target for snakemake. The cluster configuration is provided as `--cluster-config cluster.json –cluster {cluster system params}` , where at `{}` you must include the command to submit jobs you would normally use, but giving it variables that you defined at `cluster.json`. Finally, the target for snakemake should be included at the end of the command and in this case is extend_reference. Using LSF, as queue management system, we ran the following command to run this workflow:



```{bash}
snakemake --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" --use-conda -k -j 50 gffcompare/extended_ref_annotation.gtf
```

