import argparse
import csv


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--gff", required=True,
		help="Input hg38 gff file.")

	parser.add_argument(
		"--outfile", required=True,
		help="Output gff file.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	chroms = {
		"NC_000001.11": "chr1",
		"NC_000002.12": "chr2",
		"NC_000003.12": "chr3",
		"NC_000004.12": "chr4",
		"NC_000005.10": "chr5",
		"NC_000006.12": "chr6",
		"NC_000007.14": "chr7",
		"NC_000008.11": "chr8",
		"NC_000009.12": "chr9",
		"NC_000010.11": "chr10",
		"NC_000011.10": "chr11",
		"NC_000012.12": "chr12",
		"NC_000013.11": "chr13",
		"NC_000014.9": "chr14",
		"NC_000015.10": "chr15",
		"NC_000016.10": "chr16",
		"NC_000017.11": "chr17",
		"NC_000018.10": "chr18",
		"NC_000019.10": "chr19",
		"NC_000020.11": "chr20",
		"NC_000021.9": "chr21",
		"NC_000022.11": "chr22",
		"NC_000023.11": "chrX",
		"NC_000024.10": "chrY"}

	with open(args.gff, "r") as f:
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
