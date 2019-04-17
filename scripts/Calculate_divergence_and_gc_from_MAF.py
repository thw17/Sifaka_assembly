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
		"--species1", required=True,
		help="REQUIRED. Name of species/assembly (must match exactly) "
		"in the MAF sequence lines that the BED targets come from).")

	parser.add_argument(
		"--species2", required=True,
		help="REQUIRED. Name of second species/assembly (must match exactly) "
		"in the MAF sequence lines to use (note: this is NOT the species from "
		"the BED file).")

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

			max_idx = len(dict_starts[maf_record.chrom]) - 1

			while True:
				# Check if at end of targets for chromosome
				if tmp_idx > max_idx:
					break

				# Get BED coords
				tmp_bed_start = dict_starts[maf_record.chrom][tmp_idx]
				tmp_bed_stop = dict_data[(maf_record.chrom, tmp_bed_start)]

				# Case 1: Bed target start and stop are both after the MAF block stop
				# Break the loop; done with current MAF block
				# Note with BED records, subtract 1 from stop coordinate (half open)
				if maf_record.stop <= tmp_bed_start:
					break

				# Case 2: Bed target stop is before MAF block start
				# Advance the loop to get the next BED target
				elif maf_record.start >= tmp_bed_stop:
					tmp_idx += 1
					continue

				# Case 3: Targets overlap and MAF start is greater than or equal to
				# the BED start
				elif maf_record.start >= tmp_bed_start:
					# Case 3a: MAF stop less than or equal to BED stop
					# In this case, the MAF block is completely contained
					# within the BED target and can be processed as is
					if maf_record.stop <= tmp_bed_stop:
						seqs_compared = comp_seq(maf_record.seq1, maf_record.seq2)
						temp_diffs = [x for x in seqs_compared if x != "x"]
						seq1_gc = get_gc_count(maf_record.seq1)
						seq2_gc = get_gc_count(maf_record.seq2)
						# add sites and differences
						data_dict[(maf_record.chrom, tmp_bed_start)][1] += len(temp_diffs)
						data_dict[(maf_record.chrom, tmp_bed_start)][2] += sum(temp_diffs)
						# add GC
						data_dict[(maf_record.chrom, tmp_bed_start)][3] += seq1_gc
						data_dict[(maf_record.chrom, tmp_bed_start)][3] += seq2_gc
						# MAF block done; advance to next block
						break
					# Case 3b: MAF stop greater than BED stop
					# In this case, the MAF block is larger than the BED target
					# and needs to be subset (from the end, not from the beginning
					# of the sequence)
					elif maf_record.stop > tmp_bed_stop:
						# Adjust sequence (shorten), but take gaps into account
						# when indexing
						temp_idx1 = -1 * (maf_record.stop - tmp_bed_stop)
						while True:
							end_seq1 = maf_record.seq1[temp_idx1:]
							gap_count = 0
							for x in end_seq1:
								if x == "-":
									gap_count += 1
							if gap_count == 0:
								break
							else:
								temp_idx1 -= gap_count
						adj_seq1 = maf_record.seq1[:temp_idx1]
						adj_seq2 = maf_record.seq2[:temp_idx1]
						# Process adjusted sequences
						seqs_compared = comp_seq(adj_seq1, adj_seq2)
						temp_diffs = [x for x in seqs_compared if x != "x"]
						seq1_gc = get_gc_count(adj_seq1)
						seq2_gc = get_gc_count(adj_seq2)
						# add sites and differences
						data_dict[(maf_record.chrom, tmp_bed_start)][1] += len(temp_diffs)
						data_dict[(maf_record.chrom, tmp_bed_start)][2] += sum(temp_diffs)
						# add GC
						data_dict[(maf_record.chrom, tmp_bed_start)][3] += seq1_gc
						data_dict[(maf_record.chrom, tmp_bed_start)][3] += seq2_gc
						# MAF record IS NOT done; keep MAF and advance to next BED
						# target
						tmp_idx1 += 1
						continue

				# Case 4: Targets overlap and MAF start is less than BED start
				elif maf_record.start < tmp_bed_start:
					# First, adjust sequences' starts

				# Case 5: There's a condition I missed. Throw an error
				else:
					raise RuntimeError("Tim missed a condition")





if __name__ == "__main__":
	main()
