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

# macs2 input files
treatment = snakemake.input.treatment
control = snakemake.input.control

# extract params
# TODO perform parameter checks
odir = snakemake.params.get("odir", ".")
prefix = snakemake.params.get("prefix", "")
maxpeaks = snakemake.params.get("maxpeaks", 100000)
pcut = snakemake.params.get("pcut", 0.05)
genome = snakemake.params.get("genome", "hs")

# get log file
log_file = snakemake.log

# the main macs2 output file
macs_out = odir + "/" + prefix + "_peaks.narrowPeak"
# the final output file we use in the end
out_file = snakemake.output

# create additional parameters for macs2
cmd_params = "--outdir {o} -n {pr} -p {pv} -g {g}".format(g=genome, o=odir, pr=prefix, pv=pcut)
cmd = "macs2 callpeak " + cmd_params
# add treatment and control
cmd = cmd + " -t {treat} -c {ctrl} &> {log}".format(treat=treatment, ctrl=control, log=log_file)

# execut macs2 command
shell(cmd)

# a final command to create zipped regionpeak files as input for the IDR script
#cmd = "sort -k 8nr,8nr {mout} | head -n {mp} > {of}.tmp && gzip -c {of}.tmp > {of} && rm {of}.tmp".format(mout=macs_out,mp=maxpeaks,of=out_file)
#shell(cmd)

