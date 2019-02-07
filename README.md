# A snakemake pipeline for ChIP-sequencing data

## Introduction

## Workflow description

## How to get started

The following command will map, merge and filter reads for all samples as well
as downsample all bam files to 20e6 reads using a cluster computation (qsub).
It will further call peaks for all samples using the macs2 peak caller.

> NOTE: The use of '--use-conda' is mandatory if not all tools required by the 
> analysis pipeline are installed on the host system. ALso, this ensures that
> all tools which are run are used in the same version as in the paper.

To obtain only peak calls for all samples, filtered for encode blacklisted peaks,
replace the keyword 'all' at the end of the snakemake call with 'filtered_peaks_all'

```{bash} 
nohup nice snakemake --use-conda -u cluster.config --jobs=100 --local-cores=1 --cluster "qsub -pe smp {threads} -hard -l job_mem={resources.mem_mb}M \
      -q {cluster.q} -cwd -V -o {log} -e {log} -N {cluster.N}" all &
```
