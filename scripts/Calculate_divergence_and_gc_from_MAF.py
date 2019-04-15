from __future__ import print_function
from __future__ import division
import argparse
import collections
import math

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--bed", required=True,
		help="REQUIRED. Input BED file containing target regions to include.")

	parser.add_argument(
		"--maf", required=True,
		help="REQUIRED. Input MAF alignment from which to calculate divergence "
		"and GC content in targets from BED")

	parser.add_argument(
		"--output", required=True,
		help="REQUIRED. Full path to and name of output file.")

	args = parser.parse_args()
	return args


def check_nucleotide(letter):
	if letter not in "aAcCgGtT":
		return False
	else:
		return True


def comp_seq(s1, s2):
	diff_list = []
	for idx, i in enumerate(s1):
		j = s2[idx]
		if check_nucleotide(i) is False:
			diff_list.append("x")
		elif check_nucleotide(j) is False:
			diff_list.append("x")
		elif i.lower() == j.lower():
			diff_list.append(0)
		else:
			diff_list.append(1)
	return diff_list


def parsed_bed(bed_line_list):
	bed_chrom = bed_line_list[0]
	bed_start = int(bed_line_list[1])
	bed_stop = int(bed_line_list[2])
	bed_block = collections.namedtuple(
		"bed_block",
		["chrom", "start", "stop"])
	return bed_block(
		chrom=bed_chrom,
		start=bed_start,
		stop=bed_stop)


def jc69(p_dist):
	if p_dist == 0:
		return 0
	else:
		term1 = 1.0 - ((4.0 / 3.0) * p_dist)
		jc_dist = -0.75 * math.log(term1)
		return jc_dist


def main():
	args = parse_args()

	# dictionary to hold (from BED) chrom : [startsite1, startsite2, startsite3, ...]
	dict_starts = {}

	# dictionary to hold (from BED) target specific data
	# (chrom, start) : [target_len, # sites in MAF, # diffs, % GC sp1, % GC sp2]
	dict_data = {}

	# Process bed file
	with open(args.bed, "r") as f:
		for line in f:
			# parse line
			bed_record = line.strip().split()
			if bed_record == []:
				continue
			if bed_record[0][0] == "#":
				continue
			line_data = parsed_bed(bed_record)
			# add start to chrom key in dict_starts
			if line_data.chrom not in dict_starts:
				dict_starts[line_data.chrom] = [line_data.start]
			else:
				dict_starts[line_data.chrom].append(line_data.start)

			# add target to dict_data
			dict_data[(line_data.chrom, line_data.start)] = [
				line_data.stop - line_data.start, 0, 0, 0, 0]

	# Process MAF file
	with open(args.maf, "r") as maffile:
		# initialize
		maf_line_counter = 0
		maf_record = mafblock(maffile, maf_line_counter)

		# mafblock returns "Done" when it hits EOF. If this happens before anything
		# else, the file is empty.
		if maf_record == "Done":
			raise ValueError("Empty MAF file")

		# Begin looping through MAF file. First record already staged
		maf_finished = False
		while True:
			# get index of first bed start that comes before MAF block start
			tmp_idx = bisect(dict_starts[maf_record.chrom], maf_record.start)
			if tmp_idx > 0:
				tmp_idx = tmp_idx - 1

			while True:
				tmp_bed_start = dict_starts[maf_record.chrom][tmp_idx]
				tmp_bed_stop = dict_data[(maf_record.chrom, tmp_bed_start)]




if __name__ == "__main__":
	main()
