# ------------------------------------------------------------------------------
# In this file we define the 'Sample' class and some methods for easy loading
# of the sample information from the sample sheet. 
# We put this here to get the main snakemake-file more readable.
#
# @author: Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Simple class to help us manage sample information
# ------------------------------------------------------------------------------
class Sample:
	# define most important properties
        id = ""
        condition = ""
        is_paired = False
        is_multilane = False
        raw_lane1 = {}
        lane1_name = ""
        raw_lane2 = {}
        lane2_name = ""

        def __str__(self):
                return "Showing sample with id {}:\n\t{}\n\t{}\n\t{}\n\t{}\n\t{}\n\t{}".format(self.id,
                        self.condition, self.is_paired, self.raw_lane1, self.raw_lane2, self.lane1_name, self.lane2_name)


# ------------------------------------------------------------------------------
# Function to load information for all samples given the corresponding
# sample sheet
#
# Format of the sample sheet (columns) should be like the following (sep="\t"):
# 1: sample id
# 2: sample group
# 3: sample condition (e.g. diet)
# 4: sample condition 2 (e.g. timepoint)
# 5: experiment type (either PE for paired-end or SE for single-end)
# 6: raw sequencing files in the format:
# LANE:file_path_forward[,LANE:file_path_reverse|LANE2:file_path_forward,LANE2:file_path_reverse]
# ------------------------------------------------------------------------------
def load_sample_info(file):
        dict = {}
        with open(file, "r") as f:
                for line in f:
                        if(line.startswith("sample_id")):
                                continue
                        el = line.strip().split("\t")

                        s = Sample()
                        s.id = el[0]
                        s.condition = el[2]
                        s.is_paired = (el[4] == "PE")

                        raw = el[5].split("|")

                        if(len(raw)>1):
                                s.is_multilane = True
                                s.raw_lane1 = raw[0].split("=")[1].split(",")
                                s.raw_lane2 = raw[1].split("=")[1].split(",")
                                s.lane1_name = raw[0].split("=")[0]
                                s.lane2_name = raw[1].split("=")[0]
                        else:
                                s.raw_lane1 = raw[0].split("=")[1].split(",")
                                s.lane1_name = raw[0].split("=")[0]

                        dict[s.id] = s
        return(dict)

# ------------------------------------------------------------------------------
# Method to map the sample ids to the corresponding files.
# ------------------------------------------------------------------------------
def id_to_merged(wildcards):
	sample = wildcards.sample
	s      = SAMPLE_SHEET[sample]
	pe     = s.is_paired
	ml     = s.is_multilane

	if(s.is_paired):
		if(s.is_multilane):
			return(config["data"]["mapped"] + "{}_{}_{}.bam".format(sample, s.lane1_name, s.lane2_name))
		else:
			return(config["data"]["mapped"] + "{}_{}_1_2.bam".format(sample, s.lane1_name))
	elif(s.is_multilane):
		return(config["data"]["mapped"] + "{}_{}_{}.bam".format(sample, s.lane1_name, s.lane2_name))
	else:
		# single end reads, no multilane
		return(config["data"]["mapped"] + "{}.bam".format(sample))


def merged_to_pairs(wildcards):
	# get sample and lane info
	sample = wildcards.sample
	s      = SAMPLE_SHEET[sample]

	if(s.is_multilane):
		lane1  = wildcards.lane1
		lane2  = wildcards.lane2
		if(s.is_paired):
			s1     = config["data"]["mapped"] + "{}_{}_1_2.bam".format(sample, lane1)
			s2     = config["data"]["mapped"] + "{}_{}_1_2.bam".format(sample, lane2)
			pairs  = [s1,s2]
		else:
			pairs  = config["data"]["mapped"] + "{}.bam".format(sample)
	elif(s.is_paired):
		print("Not implemented yet.")
	else:
		pairs = s.raw_lane1
	return(pairs)

def pairs_to_fastq(wildcards):
        # get sample and lane info
        sample = wildcards.sample
        lane   = wildcards.lane
        s      = SAMPLE_SHEET[sample]

        if(lane == s.lane1_name):
                return(s.raw_lane1)
        if(lane == s.lane2_name):
                return(s.raw_lane2)

