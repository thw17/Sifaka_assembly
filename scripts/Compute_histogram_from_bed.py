import argparse
import csv


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--bed", required=True,
		help="Input hg38 gff file.")

	parser.add_argument(
		"--outfile", required=True,
		help="Output gff file.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()
	hist_dict = {}
	with open(args.bed, "r") as f:
		reader = csv.reader(f, delimiter="\t")
		for i in reader:
			if int(i[3]) not in hist_dict:
				hist_dict[int(i[3])] = int(i[2]) - int(i[1])
			else:
				hist_dict[int(i[3])] += int(i[2]) - int(i[1])

	key_list = hist_dict.keys()
	key_list_sorted = sorted(key_list)

	with open(args.outfile, "w") as o:
		for i in key_list_sorted:
			o.write("{}\t{}\n".format(i, hist_dict[i]))


if __name__ == "__main__":
	main()
