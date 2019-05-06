# Tim Webster, 2019, University of Utah

import argparse
import collections

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--maf", required=True,
		help="REQUIRED. Input MAF alignment.")

	parser.add_argument(
		"--species", required=True,
		help="REQUIRED. Name of species/assembly (must match exactly) "
		"in the MAF sequence lines to extract bed coordinates for.")

	parser.add_argument(
		"--print_frequency", default=1000,
		help="Print update after every X number of records processed. Default is 1000.")

	parser.add_argument(
		"--output", required=True,
		help="REQUIRED. Full path to and name of output file.")

	args = parser.parse_args()
	return args


def mafblock(file, line_num, sp1):
	"""Parses MAF blocks and returns data from one target species. Assumes """
	"""a species.chromosome format"""
	sp1_chrom = None
	sp1_start = None
	sp1_length = None
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
			# print(line)
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
			else:
				continue

	if any(x is None for x in [sp1_chrom, sp1_start, sp1_length]):
		raise RuntimeError("MAF block beginning on line {} incomplete.".format(line_num))

	maf_block = collections.namedtuple(
		"maf_block",
		["chrom", "start", "stop", "length", "counter"])

	return maf_block(
		chrom=sp1_chrom,
		start=sp1_start,
		length=sp1_length,
		stop=sp1_start + sp1_length,
		counter=line_counter)

def main():
	args = parse_args()

	# Process MAF file
	print("Starting traversal through MAF file...")
	with open(args.maf, "r") as maffile:
		with open(args.output, "w") as o:
			# initialize
			maf_line_counter = 0
			maf_record_counter = 0

			maf_record = mafblock(maffile, maf_line_counter, args.species)
			if maf_record == "Done":
				raise ValueError("Empty MAF file")
			maf_line_counter = maf_record.counter
			maf_record_counter += 1

			# Begin looping through MAF file. First record already staged
			maf_finished = False
			while maf_finished is False:
				o.write("{}\t{}\t{}\n".format(
					maf_record.chrom, maf_record.start, maf_record.stop))


				if maf_record_counter % int(args.print_frequency) == 0:
					print("{} MAF records processed...".format(maf_record_counter))

				maf_record = mafblock(maffile, maf_line_counter, args.species)
				if maf_record == "Done":
					maf_finished = True
					break
				maf_line_counter = maf_record.counter
				maf_record_counter += 1


if __name__ == "__main__":
	main()
