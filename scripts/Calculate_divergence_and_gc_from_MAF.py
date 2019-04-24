# Tim Webster, 2019, University of Utah

from __future__ import print_function
from __future__ import division
from bisect import bisect
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
		"--include_chrom_species1", nargs="*",
		help="Names of chromsomes to include from species1. Optional. All chromosmes "
		"included by default.")

	parser.add_argument(
		"--include_chrom_species2", nargs="*",
		help="Names of chromsomes to include from species2. Optional. All chromosmes "
		"included by default.")

	parser.add_argument(
		"--maf", required=True,
		help="REQUIRED. Input MAF alignment from which to calculate divergence "
		"and GC content in targets from BED")

	parser.add_argument(
		"--print_frequency", default=100,
		help="Print update after every X number of lines. Default is 100.")

	parser.add_argument(
		"--output", required=True,
		help="REQUIRED. Full path to and name of output file.")

	args = parser.parse_args()
	return args


def mafblock(file, line_num, sp1, sp2):
	"""Parses MAF blocks and returns data from two target species. Assumes """
	"""a species.chromosome format for both species"""
	sp1_chrom = None
	sp2_chrom = None
	sp1_start = None
	sp1_length = None
	sp1_seq = None
	sp2_seq = None
	line_counter = None

	skip_line = True
	line_counter = line_num - 1
	while skip_line is True:
		line = file.readline().strip().split()
		line_counter += 1
		# Check for EOF
		if not line:
			return "Done"
		# alignment blocks start with a line that begins with "a"
		if line[0] != "a":
			continue
		else:
			print(line)
			skip_line = False

	same_block = True
	while same_block is True:
		line = file.readline().strip().split()
		line_counter += 1
		if line == []:
			break
		if line[0] != "s":
			continue
		else:
			tmp_species = line[1].split(".")[0]
			tmp_scaffold = line[1].split(".")[1]
			if tmp_species == sp1:
				sp1_chrom = tmp_scaffold
				sp1_start = int(line[2])
				sp1_length = int(line[3])
				sp1_seq = line[6]
			elif tmp_species == sp2:
				sp2_chrom = tmp_scaffold
				sp2_start = int(line[2])
				sp2_length = int(line[3])
				sp2_seq = line[6]
			else:
				continue

	if any(x is None for x in [sp1_chrom, sp2_chrom, sp1_start, sp1_length, sp1_seq, sp2_seq]):
		raise RuntimeError("MAF block beginning on line {} incomplete.".format(line_num))

	maf_block = collections.namedtuple(
		"maf_block",
		["chrom", "chrom2", "start", "stop", "length", "seq1", "seq2", "counter"])

	return maf_block(
		chrom=sp1_chrom,
		chrom2=sp2_chrom,
		start=sp1_start,
		length=sp1_length,
		stop=sp1_start + sp1_length,
		seq1=sp1_seq,
		seq2=sp2_seq,
		counter=line_counter)


def check_nucleotide(letter):
	if letter not in "aAcCgGtT":
		return False
	else:
		return True


def get_gc_count(sequence):
	gc_count = 0
	total_count = 0
	for i in sequence:
		if check_nucleotide(i) is False:
			continue
		if i.lower() in ["g", "c"]:
			gc_count += 1
			total_count += 1
		else:
			total_count += 1
	return (gc_count, total_count)


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
	# (chrom, start) : [target_len, # sites in MAF, # diffs, sp1 GC count,
	# sp1 GC number of total sites, sp2 GC count, sp2 GC number of total sites]
	# use an OrderedDict for writing to output
	data_dict = collections.OrderedDict()

	# Process bed file
	print("Processing bed file...")
	with open(args.bed, "r") as f:
		for line in f:
			# parse line
			bed_record = line.strip().split()
			if bed_record == []:
				continue
			if bed_record[0][0] == "#":
				continue
			line_data = parsed_bed(bed_record)
			# check if chrom is included
			if args.include_chrom_species1 is not None:
				if line_data.chrom not in args.include_chrom_species1:
					continue
			# add start to chrom key in dict_starts
			if line_data.chrom not in dict_starts:
				dict_starts[line_data.chrom] = [line_data.start]
			else:
				dict_starts[line_data.chrom].append(line_data.start)

			# add target to data_dict
			data_dict[(line_data.chrom, line_data.start)] = [
				line_data.stop - line_data.start, 0, 0, 0, 0, 0, 0]
	print("Processing bed file complete.")

	# Process MAF file
	print("Starting traversal through MAF file...")
	with open(args.maf, "r") as maffile:
		# initialize
		maf_line_counter = 0
		maf_record_counter = 0
		# maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
		# maf_record_counter += 1

		if args.include_chrom_species1 is None and args.include_chrom_species2 is None:
			maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
			if maf_record == "Done":
				raise ValueError("Empty MAF file")
			maf_line_counter = maf_record.counter
			maf_record_counter += 1
		elif args.include_chrom_species1 is not None:
			skip_chrom = True
			if args.include_chrom_species2 is None:
				while skip_chrom is True:
					maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
					if maf_record == "Done":
						raise ValueError("Either empty MAF file or included chromosomes not in file")
					maf_line_counter = maf_record.counter
					maf_record_counter += 1
					if maf_record.chrom not in args.include_chrom_species1:
						continue
					else:
						skip_chrom = False
			else:
				while skip_chrom is True:
					maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
					if maf_record == "Done":
						raise ValueError("Either empty MAF file or included chromosomes not in file")
					maf_line_counter = maf_record.counter
					maf_record_counter += 1
					if maf_record.chrom not in args.include_chrom_species1:
						continue
					else:
						if maf_record.chrom2 not in args.include_chrom_species2:
							continue
						else:
							skip_chrom = False
		elif args.include_chrom_species2 is not None:
			skip_chrom = True
			while skip_chrom is True:
				maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
				if maf_record == "Done":
					raise ValueError("Either empty MAF file or included chromosomes not in file")
				maf_line_counter = maf_record.counter
				maf_record_counter += 1
				if maf_record.chrom2 not in args.include_chrom_species2:
					continue
				else:
					skip_chrom = False

		# mafblock returns "Done" when it hits EOF. If this happens before anything
		# else, the file is empty.
		if maf_record == "Done":
			raise ValueError("Empty MAF file")

		maf_line_counter = maf_record.counter

		# Begin looping through MAF file. First record already staged
		maf_finished = False
		while maf_finished is False:
			# get index of first bed start that comes before MAF block start
			if maf_record.chrom in dict_starts:
				advance_maf = False
			else:
				advance_maf = True

			if advance_maf is False:
				bed_idx = bisect(dict_starts[maf_record.chrom], maf_record.start)
				if bed_idx > 0:
					bed_idx = bed_idx - 1

				max_idx = len(dict_starts[maf_record.chrom]) - 1

				while True:
					# Check if at end of targets for chromosome
					if bed_idx > max_idx:
						break

					# Get BED coords
					tmp_bed_start = dict_starts[maf_record.chrom][bed_idx]
					tmp_bed_stop = data_dict[(maf_record.chrom, tmp_bed_start)][0] + tmp_bed_start

					# Case 1: Bed target start and stop are both after the MAF block stop
					# Break the loop; done with current MAF block
					# Note with BED records, subtract 1 from stop coordinate (half open)
					if maf_record.stop <= tmp_bed_start:
						print("case1")
						break

					# Case 2: Bed target stop is before MAF block start
					# Advance the loop to get the next BED target
					elif maf_record.start >= tmp_bed_stop:
						print("case2")
						bed_idx += 1
						continue

					# Case 3: Targets overlap and MAF start is greater than or equal to
					# the BED start
					elif maf_record.start >= tmp_bed_start:
						# Case 3a: MAF stop less than or equal to BED stop
						# In this case, the MAF block is completely contained
						# within the BED target and can be processed as is
						if maf_record.stop <= tmp_bed_stop:
							print("case3a")
							seqs_compared = comp_seq(maf_record.seq1, maf_record.seq2)
							temp_diffs = [x for x in seqs_compared if x != "x"]
							seq1_gc = get_gc_count(maf_record.seq1)
							seq2_gc = get_gc_count(maf_record.seq2)
							# add sites and differences
							data_dict[(maf_record.chrom, tmp_bed_start)][1] += len(temp_diffs)
							data_dict[(maf_record.chrom, tmp_bed_start)][2] += sum(temp_diffs)
							# add GC
							data_dict[(maf_record.chrom, tmp_bed_start)][3] += seq1_gc[0]
							data_dict[(maf_record.chrom, tmp_bed_start)][4] += seq1_gc[1]
							data_dict[(maf_record.chrom, tmp_bed_start)][5] += seq2_gc[0]
							data_dict[(maf_record.chrom, tmp_bed_start)][6] += seq2_gc[1]
							# MAF block done; advance to next block
							break
						# Case 3b: MAF stop greater than BED stop
						# In this case, the MAF block is larger than the BED target
						# and needs to be subset (from the end, not from the beginning
						# of the sequence)
						elif maf_record.stop > tmp_bed_stop:
							print("case3b")
							# Adjust sequence (shorten), but take gaps into account
							# when indexing
							# Here, temp_idx1 will always be used for negative indexing
							# (i.e., counting from the end)
							temp_idx1 = -1 * (maf_record.stop - tmp_bed_stop)
							orig_idx = temp_idx1
							while True:
								end_seq1 = maf_record.seq1[temp_idx1:]
								gap_count = 0
								for x in end_seq1:
									if x == "-":
										gap_count += 1
								if len(end_seq1) - gap_count == -1 * (orig_idx):
									break
								# if gap_count == 0:
								# 	break
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
							data_dict[(maf_record.chrom, tmp_bed_start)][3] += seq1_gc[0]
							data_dict[(maf_record.chrom, tmp_bed_start)][4] += seq1_gc[1]
							data_dict[(maf_record.chrom, tmp_bed_start)][5] += seq2_gc[0]
							data_dict[(maf_record.chrom, tmp_bed_start)][6] += seq2_gc[1]
							# MAF record IS NOT done; keep MAF and advance to next BED
							# target
							bed_idx += 1
							continue

					# Case 4: Targets overlap and MAF start is less than BED start
					elif maf_record.start < tmp_bed_start:
						# First, adjust sequences' starts while taking into account gaps
						# - need to bring MAF start for seq1 up to the beginning of the bed region
						# temp_idx2 will be used for indexing from the start (positive numbers only)
						temp_idx2 = tmp_bed_start - maf_record.start
						orig_idx = temp_idx2
						while True:
							beg_seq1 = maf_record.seq1[:temp_idx2]
							gap_count = 0
							for x in beg_seq1:
								if x == "-":
									gap_count +=1
							if temp_idx2 - gap_count == orig_idx:
								break
							else:
								temp_idx2 += gap_count
						adj_seq1 = maf_record.seq1[:temp_idx2]
						adj_seq2 = maf_record.seq2[:temp_idx2]
						# Case 4a: MAF block end is less than or equal to BED end
						# No further adjustments are necessary for MAF block, calculate stats
						if maf_record.stop <= tmp_bed_stop:
							print("case4a")
							# Process adjusted sequences
							seqs_compared = comp_seq(adj_seq1, adj_seq2)
							temp_diffs = [x for x in seqs_compared if x != "x"]
							seq1_gc = get_gc_count(adj_seq1)
							seq2_gc = get_gc_count(adj_seq2)
							# add sites and differences
							data_dict[(maf_record.chrom, tmp_bed_start)][1] += len(temp_diffs)
							data_dict[(maf_record.chrom, tmp_bed_start)][2] += sum(temp_diffs)
							# add GC
							data_dict[(maf_record.chrom, tmp_bed_start)][3] += seq1_gc[0]
							data_dict[(maf_record.chrom, tmp_bed_start)][4] += seq1_gc[1]
							data_dict[(maf_record.chrom, tmp_bed_start)][5] += seq2_gc[0]
							data_dict[(maf_record.chrom, tmp_bed_start)][6] += seq2_gc[1]
							# MAF record is done; move to next MAF
							break
						# Case 4b: MAF block end is greater than BED; subset MAF before
						# calculating stats and then advance BED target
						elif maf_record.stop > tmp_bed_stop:
							print("case4b")
							temp_idx1 = -1 * (maf_record.stop - tmp_bed_stop)
							orig_idx = -1 * temp_idx1
							while True:
								end_seq1 = maf_record.seq1[temp_idx1:]
								gap_count = 0
								for x in end_seq1:
									if x == "-":
										gap_count += 1
								if orig_idx + gap_count == len(end_seq1):
									break
								else:
									temp_idx1 -= gap_count
							# Subset from the *already adjusted sequences*
							adj_seq1 = adj_seq1[:temp_idx1]
							adj_seq2 = adj_seq2[:temp_idx1]
							# Process adjusted sequences
							seqs_compared = comp_seq(adj_seq1, adj_seq2)
							temp_diffs = [x for x in seqs_compared if x != "x"]
							seq1_gc = get_gc_count(adj_seq1)
							seq2_gc = get_gc_count(adj_seq2)
							# add sites and differences
							data_dict[(maf_record.chrom, tmp_bed_start)][1] += len(temp_diffs)
							data_dict[(maf_record.chrom, tmp_bed_start)][2] += sum(temp_diffs)
							# add GC
							data_dict[(maf_record.chrom, tmp_bed_start)][3] += seq1_gc[0]
							data_dict[(maf_record.chrom, tmp_bed_start)][4] += seq1_gc[1]
							data_dict[(maf_record.chrom, tmp_bed_start)][5] += seq2_gc[0]
							data_dict[(maf_record.chrom, tmp_bed_start)][6] += seq2_gc[1]
							# MAF record IS NOT done; keep MAF and advance to next BED
							# target
							bed_idx += 1
							continue

					# Case 5: There's a condition I missed. Throw an error
					else:
						print("case5")
						raise RuntimeError("Tim missed a condition")

			if maf_record_counter % int(args.print_frequency) == 0:
				print("{} MAF records processed...".format(maf_record_counter))
			# Advance MAF record and end if done
			# case of considering all chromosomes
			if args.include_chrom_species1 is None and args.include_chrom_species2 is None:
				maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
				if maf_record == "Done":
					maf_finished = True
					break
				maf_line_counter = maf_record.counter
				maf_record_counter += 1
			elif args.include_chrom_species1 is not None:
				skip_chrom = True
				if args.include_chrom_species2 is None:
					while skip_chrom is True:
						maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
						if maf_record == "Done":
							maf_finished = True
							break
						maf_line_counter = maf_record.counter
						maf_record_counter += 1
						if maf_record.chrom not in args.include_chrom_species1:
							continue
						else:
							skip_chrom = False
				else:
					while skip_chrom is True:
						maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
						if maf_record == "Done":
							maf_finished = True
							break
						maf_line_counter = maf_record.counter
						maf_record_counter += 1
						if maf_record.chrom not in args.include_chrom_species1:
							continue
						else:
							if maf_record.chrom2 not in args.include_chrom_species2:
								continue
							else:
								skip_chrom = False
			elif args.include_chrom_species2 is not None:
				skip_chrom = True
				while skip_chrom is True:
					maf_record = mafblock(maffile, maf_line_counter, args.species1, args.species2)
					if maf_record == "Done":
						maf_finished = True
						break
					maf_line_counter = maf_record.counter
					maf_record_counter += 1
					if maf_record.chrom2 not in args.include_chrom_species2:
						continue
					else:
						skip_chrom = False

	# Output results
	with open(args.output, "w") as f:
		f.write("chrom\tstart\tlength\tnum_sites\tdiffs\tp_dist\tJC69\tfrac_gc_sp1\tfrac_gc_sp2\n")
		for key in data_dict:
			tmp_list = data_dict[key]
			if tmp_list[1] == 0:
				p_dist = 0
				jc = 0
				frac_gc_sp1 = 0
				frac_gc_sp2 = 0
			else:
				p_dist = tmp_list[2] / tmp_list[1]
				jc = jc69(tmp_list[2] / tmp_list[1])
				frac_gc_sp1 = tmp_list[3] / tmp_list[4]
				frac_gc_sp2 = tmp_list[5] / tmp_list[6]
			f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
				key[0],
				key[1],
				tmp_list[0],
				tmp_list[1],
				tmp_list[2],
				p_dist,
				jc,
				frac_gc_sp1,
				frac_gc_sp2))


if __name__ == "__main__":
	main()
