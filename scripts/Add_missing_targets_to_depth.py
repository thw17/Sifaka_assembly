import argparse
import collections

def parse_args():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--input_depth", required=True,
		help="REQUIRED. Input per target BED file with depth stats.")

	parser.add_argument(
		"--input_maf_stats", required=True,
		help="REQUIRED. Input per target BED file with maf stats.")

	parser.add_argument(
		"--output", required=True,
		help="REQUIRED. Output BED file.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	depth_targets = set()
	with open(args.input_depth, "r") as f:
		for line in f:
			line1 = line.strip().split()
			if line1 == []:
				continue
			if line1[0] == "chrom":
				continue
			depth_targets.add((line1[0], line1[1]))

	missing_targets = collections.OrderedDict()
	with open(args.input_maf_stats, "r") as f:
		for line in f:
			line1 = line.strip().split()
			if line1 == []:
				continue
			if line1[0] == "chrom":
				continue
			chrom_start = (line1[0], line1[1])
			if chrom_start not in depth_targets:
				missing_targets[chrom_start] = [line1[0], line1[1], int(line1[1]) + int(line1[2]), line1[2], 0]

	with open(args.input_depth, "r") as f:
		with open(args.output, "w") as o:
			for line in f:
				line1 = line.strip().split()
				if line1 == []:
					continue
				out_line = "\t".join(line1)
				o.write("{}\n".format(out_line))
			for item in missing_targets:
				out_line = "\t".join([str(x) for x in missing_targets[item]])
				o.write("{}\n".format(out_line))


if __name__ == "__main__":
	main()
