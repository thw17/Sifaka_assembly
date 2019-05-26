import argparse

def parse_args():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--input", required=True,
		help="REQUIRED. Input per target BED file with maf stats.")

	parser.add_argument(
		"--output", required=True,
		help="REQUIRED. Output BED file.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	with open(args.input, "r") as f:
		with open(args.output, "w") as o:
			for line in f:
				line1 = line.strip().split()
				if line1 == []:
					continue
				if line1[0] == "chrom":
					line1.append("chrom_start")
					out_line = "\t".join(line1)
					o.write("{}\n".format(out_line))
					continue
				line1.append("{}_{}".format(line1[0], line1[1]))
				out_line = "\t".join(line1)
				o.write("{}\n".format(out_line))


if __name__ == "__main__":
	main()
