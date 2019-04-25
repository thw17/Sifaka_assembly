import argparse
import csv


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--bed", required=True,
		help="Input Mmul8.0.1 bed file.")

	parser.add_argument(
		"--outfile", required=True,
		help="Output bed file.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	chroms = {
		"NC_027893.1": "chr1",
		"NC_027894.1": "chr2",
		"NC_027895.1": "chr3",
		"NC_027896.1": "chr4",
		"NC_027897.1": "chr5",
		"NC_027898.1": "chr6",
		"NC_027899.1": "chr7",
		"NC_027900.1": "chr8",
		"NC_027901.1": "chr9",
		"NC_027902.1": "chr10",
		"NC_027903.1": "chr11",
		"NC_027904.1": "chr12",
		"NC_027905.1": "chr13",
		"NC_027906.1": "chr14",
		"NC_027907.1": "chr15",
		"NC_027908.1": "chr16",
		"NC_027909.1": "chr17",
		"NC_027910.1": "chr18",
		"NC_027911.1": "chr19",
		"NC_027912.1": "chr20",
		"NC_027913.1": "chrX",
		"NC_027914.1": "chrY"}

	with open(args.bed, "r") as f:
		with open(args.outfile, "w") as o:
			reader = csv.reader(f, delimiter="\t")
			writer = csv.writer(o, delimiter="\t")
			for i in reader:
				temp_line = i
				if i[0] in chroms:
					temp_line[0] = chroms[i[0]]
				writer.writerow(temp_line)


if __name__ == "__main__":
	main()
