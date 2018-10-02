__author__ = "Johann Hawe"
__copyright__ = "Copyright 2018, Johann Hawe"
__email__ = "johann.hawe@helmholtz-muenchen.de"
#__license__ = "MIT"


from snakemake.shell import shell


shell("samtools view -b -s0.5 {snakemake.params} -@{snakemake.threads} -U {snakemake.output[0]} {snakemake.input[0]} > {snakemake.output[1]} && samtools index {snakemake.output[0]} && samtools index {snakemake.output[1]}")
