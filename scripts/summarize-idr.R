# define cutoff
cutoff <- -log10(snakemake@params$pcut)

# get the output file
out_file <- snakemake@output[[1]]

# get sample name
sample <- snakemake@params$sample

# get the result files and read in the data
fselfcons_rep1 <- snakemake@input[["rep1_selfconsistency"]]
fselfcons_rep2 <- snakemake@input[["rep2_selfconsistency"]]
fmerge_selfcons <- snakemake@input[["merged_selfconsistency"]]
frep1_rep2 <- snakemake@input[["rep1_rep2"]]

files <- c(fselfcons_rep1, fselfcons_rep2, fmerge_selfcons, frep1_rep2)
names(files) <- c("rep1_selfconsistency", "rep2_selfconsistency", 
			"merged_selfconsistency", "rep1_rep2")

# in the end we want to take the global FDR which is in the 12th col
# and get the number of peaks meating the specified cutoff
npeaks <- sapply(files, function(f) {
	peaks <- read.table(f, header=T, sep="\t")
	npeaks <- sum(peaks[,12] > cutoff)
	npeaks
})
names(npeaks) <- names(files)

# now calculate the IDR measure indicating (non)-reproducible replicates
r1 <- npeaks["rep1_selfconsistency"] / npeaks["rep2_selfconsistency"]
r2 <- npeaks["merged_selfconsistency"] / npeaks["rep1_rep2"]

# if both r1 and r2 are > 2/<0.5 --> low reproducibility
is_low_repr <- (r1 > 2 | r1 < 0.5) & (r2 > 2 | r2 < 0.5)
passed <- !is_low_repr

# now report the measures and write them to a file
out <- rbind(c(sample, npeaks, r1, r2, passed))
colnames(out) <- c("sample", names(npeaks), "n1_n2","np_nt", "passed")

write.table(file=out_file, sep="\t", col.names=T, row.names=F, quote=F, 
		out)
