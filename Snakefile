configfile: "./config.json"

# lists containing used variables for generating file names
# suffixes for replicates
RSU = [config["replicates"]["suffix1"], config["replicates"]["suffix2"]]
# suffixes for split bam files in idr analysis
SSU = [config["split_bam"]["suffix1"],config["split_bam"]["suffix2"]]

include: "sm_includes/process_bam.sm"
include: "sm_includes/call_peaks.sm"
#include: "sm_includes/idr.sm"
