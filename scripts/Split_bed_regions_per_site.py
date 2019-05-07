# Tim Webster, 2019, University of Utah

import argparse
import csv


def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--bed", required=True,
		help="REQUIRED. Input bedfile.")

	parser.add_argument(
		"--output", required=True,
		help="REQUIRED. Full path to and name of output file.")

	args = parser.parse_args()
	return args


def main():
	args = parse_args()

	with open(args.bed, "r") as f:
		with open(args.output, "w") as o:
			for i in csv.reader(f, delimiter="\t"):
				if i[0][0] != "#":
					start = int(i[1])
					stop = int(i[2])
					if stop - start == 1:
						o.write(
							"{}\n".format(
								"\t".join(i)))
					else:
						stop_2 = start + 1
						start_2 = start
						while stop_2 <= stop:
							output_line = [i[0], start_2, stop_2] + i[3:]
							output_line = [str(x) for x in output_line]
							o.write(
								"{}\n".format(
									"\t".join(output_line)))
							stop_2 += 1
							start_2 += 1


if __name__ == "__main__":
	main()
