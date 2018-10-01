# A snakemake pipeline for ChIP-sequencing data

## Introduction

## Workflow description

## How to get started

The following command will map, merge and filter reads for all samples as well as downsample all bam files to 20e6 reads using a cluster computation (qsub):
nohup nice snakemake -u cluster.config --jobs=100 --local-cores=1 --cluster "qsub -pe smp {threads} -hard -l job_mem={resources.mem_mb}M \
      -q {cluster.q} -cwd -V -o {log} -e {log} -N {cluster.N}" process_all &
