# -------------------------------------------------------------------------------
#' Remove the blacklisted encode peaks from the called peaks of a sample.
#' Removes all peaks overlapping with a blacklisted peak.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# -------------------------------------------------------------------------------
sink(file=snakemake@log[[1]], append=F, type="output", split=T)

# -------------------------------------------------------------------------------
print("Loading libraries and scripts.")
# -------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

# -------------------------------------------------------------------------------
print("Getting snakemake params.")
# -------------------------------------------------------------------------------
# inputs
fpeaks <- snakemake@input$peaks
fblacklist <- snakemake@input$blacklist

# output
fout <- snakemake@output[[1]]

# -------------------------------------------------------------------------------
print("Processing peaks.")
# -------------------------------------------------------------------------------
peaks <- read.table(fpeaks, header=F, sep="\t")
peaks <- with(peaks, GRanges(V1 ,IRanges(V2, V3)))
blacklist <- import(fblacklist)

# get indices of overlapping peaks
idxs <- overlapsAny(peaks, blacklist)
peaks <- peaks[idxs]

print(paste0("Removed ", sum(!idxs), " peaks from peak list."))

# export the new peak list
# export(peaks, fout) has bug for empty peak lists
write.table(cbind(as.character(seqnames(peaks)), start(peaks), end(peaks)),
            col.names=F, row.names=F, quote=F, sep="\t", file=fout)

# -------------------------------------------------------------------------------
print("SessionInfo:")
# -------------------------------------------------------------------------------
sessionInfo()

sink()

