configfile: "./config.json"

# lists containing used variables for generating file names
# suffixes for replicates
RSU = [config["replicates"]["suffix1"], config["replicates"]["suffix2"]]
# suffixes for split bam files in idr analysis
SSU = [config["split_bam"]["suffix1"],config["split_bam"]["suffix2"]]


# ------------------------------------------------------------------------------
# include all rules from the sub-processes and the methods file
# ------------------------------------------------------------------------------
include: "sm_includes/methods.sm"
include: "sm_includes/process_bam.sm"
include: "sm_includes/call_peaks.sm"
include: "sm_includes/mapping.sm"
#include: "sm_includes/idr.sm"

# load the samplesheet
SAMPLE_SHEET = load_sample_info(config["data"]["sample_sheet"])
print("Loaded samples:")
print(SAMPLE_SHEET.keys())

# here we define the actual samples to be processed
# TODO: make this more handy, e.g. generate list etc...
SAMPLES = ["AL11", "AL22"]#,"","",""]


# rules to exclude for cluster processing
localrules:
	process_all

# define global rules

rule process_all:
	input:
		expand(config["results"]["processed"] + "{sample}_20000000.bam", sample=SAMPLES)
