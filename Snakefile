configfile: "./configs/workflow.json"

# ------------------------------------------------------------------------------
# include all rules from the sub-processes and the methods file
# ------------------------------------------------------------------------------
include: "sm_includes/methods.py"
include: "sm_includes/process_bam.sm"
include: "sm_includes/call_peaks.sm"
include: "sm_includes/mapping.sm"
include: "sm_includes/idr.sm"
include: "sm_includes/peak_universe.sm"

# load the samplesheet
SAMPLE_SHEET = load_sample_info(config["data"]["sample_sheet"])

# here we define the actual samples to be processed
# on the long run, we want to extract this from the 
# sample sheet.
#SAMPLES = SAMPLE_SHEET.keys()
#print(SAMPLES)
# For now we just define two samples for testing:
SAMPLES = ["BH4-1", "BH5-1"]
# needed to perform IDR, TODO
SAMPLES_NOREP = ["BH4", "BH5"]

# ------------------------------------------------------------------------------
# Rules to exclude for cluster processing
# ------------------------------------------------------------------------------
localrules:
	process_all, call_peaks_all, link_replicate, 
        filtered_peaks_all, remove_blacklisted_peaks

# ------------------------------------------------------------------------------
# Define global target rules for different steps of the pipeline
# ------------------------------------------------------------------------------
rule process_all:
	input:
		expand(config["results"]["processed"] + 
                "{{sample}}_{}.bam".format(config["downsample"]["size"]), 
                sample=SAMPLES)

rule call_peaks_all:
	input:
		expand(config["results"]["peaks"] + "{sample}_peaks.narrowPeak", sample=SAMPLES)

rule filtered_peaks_all:
	input:
		expand(config["results"]["peaks"] + "{sample}_blacklistremoved.bed", sample=SAMPLES)

rule idr_all:
        input:  
                 expand(config["results"]["idr"] + "{sample}_idr_summary.txt", sample=SAMPLES_NOREP)
        output:
                config["results"]["idr"] + "idr_summary.txt"
        shell:
                "cat {input} | sort -r | uniq > {output}"

# Create summary plots for the IDR analysis.
rule plot_idr:
        input: config["results"]["idr"] + "idr_summary.txt"
        output: config["results"]["idr"] + "idr_summary.html"
        script:
                "scripts/plot-idr.Rmd"

# calls the complete pipeline (peak calling, idr analysis, still TODO)
rule all:
	input:
		config["results"]["idr"] + "idr_summary.html",
		expand(config["results"]["peaks"] + "{sample}_peaks.narrowPeak", sample=SAMPLES)
