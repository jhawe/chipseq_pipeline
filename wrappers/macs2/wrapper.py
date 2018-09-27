__author__ = "Johann Hawe"
__copyright__ = "Copyright 2018, Johann Hawe"
__email__ = "johann.hawe@helmholtz-muenchen.de"
#__license__ = "MIT"


from snakemake.shell import shell

# call peaks using the macs2 peak caller

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
maxpeaks = snakemake.params.get("maxpeaks", -1)
pcut = snakemake.params.get("pcut", -1)
qcut = snakemake.params.get("qcut", -1)

cutoff_param = " -p "
cutoff = -1

if(pcut != -1):
	cutoff = pcut
else if(pcut == -1 && qcut == -1):
	cutoff = 0.05
else:
	cutoff = qcut
	cutoff_param = " -q "

genome = snakemake.params.get("genome", "hs")

# get log file
log_file = snakemake.log

# the main macs2 output file
macs_out = odir + "/" + prefix + "_peaks.narrowPeak"
# the final output file we use in the end
out_file = snakemake.output

# create additional parameters for macs2
cmd_params = "--outdir {o} -n {pr} -g {g}".format(g=genome, o=odir, pr=prefix)
cmd_params = cmd_params + cutoff_param + "{v}".format(v=cutoff)
cmd = "macs2 callpeak " + cmd_params
# add treatment and control
cmd = cmd + " -t {treat} -c {ctrl} &> {log}".format(treat=treatment, ctrl=control, log=log_file)

# execut macs2 command
shell(cmd)
