__author__ = "Johann Hawe"
__copyright__ = "Copyright 2018, Johann Hawe"
__email__ = "johann.hawe@helmholtz-muenchen.de"
#__license__ = "MIT"

from snakemake.shell import shell

# call peaks with low stringency (params via snakemake rule)
# using macs2.

# output definition for macs2 defined in params
# needs the following values at least:
# -odir = "<macs-out-dir>"
# -prefix = "<macs-out-prefix>"

# idr input files
rep1 = snakemake.input.rep1
rep2 = snakemake.input.rep2

# extract params
# TODO perform parameter checks
#odir = snakemake.params.get("odir", ".")
#prefix = snakemake.params.get("prefix", "")
#maxpeaks = snakemake.params.get("maxpeaks", 100000)
#pcut = snakemake.params.get("pcut", 0.05)
#genome = snakemake.params.get("genome", "hs")

# get log file
log_file = snakemake.log

# the final output file
out_file = snakemake.output

# create additional parameters for idr analysis
cmd_params = "--verbose --rank {ra} -o {o} -l {l} --input-file-type {t}".format(t="narrowPeak", o=out_file, l=log_file, ra="p.value")
cmd = "idr {pars} --plot -s {r1} {r2}".format(r1=rep1, r2=rep2, pars=cmd_params)

# execut idr command
shell(cmd)
