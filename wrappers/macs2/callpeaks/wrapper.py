__author__ = "Johann Hawe"
__copyright__ = "Copyright 2018, Johann Hawe"
__email__ = "johann.hawe@helmholtz-muenchen.de"
#__license__ = "MIT"

# ------------------------------------------------------------------------------
# Call peaks using the macs2 peak caller
# ------------------------------------------------------------------------------
import sys
from sm_includes import methods
from snakemake.shell import shell

# output definition for macs2 defined in params
# needs the following values at least:
# -odir = "<macs-out-dir>"
# -prefix = "<macs-out-prefix>"

# get the sample information
# we need this to get the paired end information
SAMPLE_SHEET = methods.load_sample_info(snakemake.config["data"]["sample_sheet"])
sample = SAMPLE_SHEET[snakemake.wildcards.sample]

# macs2 input files
treatment = snakemake.input.treatment
control = snakemake.input.control

# extract params
odir = snakemake.params.get("odir", ".")
prefix = snakemake.params.get("prefix", "")
maxpeaks = snakemake.params.get("maxpeaks", -1)
pcut = snakemake.params.get("pcut", -1)
qcut = snakemake.params.get("qcut", -1)
paired = sample.is_paired

# define the cutoff param to be used (either -p or -q)
cutoff_param = " -p "
cutoff = -1

# p-value cutoff provided -> use this
if(pcut != -1):
	cutoff = pcut
# not cut-off provided, default is to use 0.05 pval cutoff
elif(pcut == -1 and qcut == -1):
	cutoff = 0.05
else:
	# use the q-value cutoff provided
	cutoff = qcut
	cutoff_param = " -q "

genome = snakemake.params.get("genome", "")
if(genome == ""):
	print("Genome has to be specified to call peaks!")
	sys.exit()

# get log file
log_file = snakemake.log

# the main macs2 output file
macs_out = odir + "/" + prefix + "_peaks.narrowPeak"
# the final output file we use in the end
out_file = snakemake.output

# create additional parameters for macs2
cmd_params = "--outdir {o} -n {pr} -g {g}".format(g=genome, o=odir, pr=prefix)
cmd_params = cmd_params + cutoff_param + "{v}".format(v=cutoff)
# explicitely tell if we have a paired-end sample!
if(paired):
	cmd_params = cmd_params + " -f BAMPE"
cmd = "macs2 callpeak " + cmd_params

# add treatment and control
cmd = cmd + " -t {treat} -c {ctrl} &> {log}".format(treat=treatment, ctrl=control, log=log_file)

# execut macs2 command
shell(cmd)
