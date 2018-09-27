__author__ = "Johann Hawe"
__copyright__ = "Copyright 2018, Johann Hawe"
__email__ = "johann.hawe@helmholtz-muenchen.de"
#__license__ = "MIT"


from snakemake.shell import shell
import datetime

# log some information
with open(snakemake.log[0], "w") as logger:
	nreads=0
	logger.write("Started at {d}\n".format(d=datetime.datetime.now()))
	logger.write("Getting number of reads from bam file.\n")

	#get total number of reads
	for line in shell("samtools view -c {snakemake.input}", iterable=True):
		nreads = line

	percent = int(snakemake.wildcards.samplesize) / int(nreads)

	msg = "Total number of reads: {r}\n".format(r=nreads)
	logger.write(msg)
	msg = "Number of reads to keep: {n}\n".format(n=snakemake.wildcards.samplesize)
	logger.write(msg)
	msg = "Percent of reads to keep: {p:.2f}\n".format(p=percent)
	logger.write(msg)

	# check whether we have > 100% -> set to 1
	if(percent > 1):
		percent = 1

	# prepare command to downsample and index files
	cmd = "samtools view -b -@{t} -s{p:.2f} {{snakemake.input}} > {{snakemake.output[0]}}".format(t=snakemake.threads, p=percent)
	cmd = cmd + " && samtools index {snakemake.output[0]}"
	shell(cmd)

	# log some more information before leaving
	for line in shell("samtools view -c {snakemake.output[0]}", iterable=True):
		logger.write("Number of reads after downsampling: {r}\n".format(r=line))
		percent = int(line) / int(nreads)
		logger.write("Percentage of reads after downsampling: {p:.2f}\n".format(p=percent))
	logger.write("Done at {d}\n".format(d=datetime.datetime.now()))

