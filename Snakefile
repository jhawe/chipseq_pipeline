configfile: "./config.json"

# ------------------------------------------------------------------------------
# include all rules from the sub-processes and the methods file
# ------------------------------------------------------------------------------
include: "sm_includes/methods.py"
include: "sm_includes/process_bam.sm"
include: "sm_includes/call_peaks.sm"
include: "sm_includes/mapping.sm"
include: "sm_includes/idr.sm"

# load the samplesheet
SAMPLE_SHEET = load_sample_info(config["data"]["sample_sheet"])

# here we define the actual samples to be processed
# TODO: make this more handy, e.g. generate list etc...
#SAMPLES = ["BH41", "BH51", "BH11", "BH12", "BH21", "BH22", "BH31", "BH32",
#"BH42","BH52","BH61","BH63","AL12","AL31","AL41","BL11","BL21","BL22-2", "BL31", "BL46", "BL51", "BL52", "BL61","BL62"]

# for testing
SAMPLES = ["AL11", "AL22"]
SAMPLES_NOREP = ["AL1", "AL2"]

# ------------------------------------------------------------------------------
# Rules to exclude for cluster processing
# ------------------------------------------------------------------------------
localrules:
	process_all, call_peaks_all

# ------------------------------------------------------------------------------
# Define global rules
# ------------------------------------------------------------------------------

rule process_all:
	input:
		expand(config["results"]["processed"] + "{{sample}}_{}.bam".format(config["downsample"]["size"]), sample=SAMPLES)

rule call_peaks_all:
	input:
		expand(config["results"]["peaks"] + "{sample}_peaks.narrowPeak", sample=SAMPLES)

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

# calls the complete pipeline (peak calling, idr analysis, ...TODO)
rule all:
	input:
		config["results"]["idr"] + "idr_summary.html",
		expand(config["results"]["peaks"] + "{sample}_peaks.narrowPeak", sample=SAMPLES)
