import argparse

def parse_args():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--input", required=True,
		help="REQUIRED. Input per target BED file with depth.")

	parser.add_argument(
		"--id", required=True,
		help="REQUIRED. ID to add to each row in same column (i.e., sample ID).")

	parser.add_argument(
		"--output", required=True,
		help="REQUIRED. Output BED file.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	with open(args.input, "r") as f:
		with open(args.output, "w") as o:
			o.write("chrom\tstart\tstop\tlength\tmean_depth\tchrom_start\tid\n")
			for line in f:
				line1 = line.strip().split()
				if line1 == []:
					continue
				line1.append("{}_{}".format(line1[0], line1[1]))
				line1.append("{}".format(args.id))
				out_line = "\t".join(line1)
				o.write("{}\n".format(out_line))


if __name__ == "__main__":
	main()
