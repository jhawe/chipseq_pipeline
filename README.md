# A snakemake pipeline for ChIP-sequencing data

## Introduction

## Workflow description

## How to get started

The following command will map, merge and filter reads for all samples as well
as downsample all bam files to 20e6 reads using a cluster computation (qsub).
It will further call peaks for all samples using the macs2 peak caller:
nohup nice snakemake --use-conda -u cluster.config --jobs=100 --local-cores=1 --cluster "qsub -pe smp {threads} -hard -l job_mem={resources.mem_mb}M \
      -q {cluster.q} -cwd -V -o {log} -e {log} -N {cluster.N}" call_peaks_all &
