import os

def main():
	bin_genes = []

	for folder_name in os.listdir("../R/data/busco_genes/"):
		for genes_faa in os.listdir("../R/data/busco_genes/" + folder_name):
			sample_bin_gene_entry = {"SampleName": folder_name[6::], "Bin": genes_faa[18:21]}

			file = open("../R/data/busco_genes/" + folder_name + "/" + genes_faa, "r")
			genes = 0

			for line in file.readlines():
				if line[0] == ">":
					genes += 1

			sample_bin_gene_entry["Genes"] = str(genes)

			bin_genes.append(sample_bin_gene_entry)

	colnames = ["SampleName", "Bin", "Genes"]

	bin_genes_tsv = ""

	for i in range(len(colnames)):
		bin_genes_tsv += colnames[i] + ("\t" if i < len(colnames) - 1 else "")

	for data_entry in bin_genes:
		bin_genes_tsv += "\n"
		for i in range(len(colnames)):
			bin_genes_tsv += str(data_entry[colnames[i]]) + ("\t" if i < len(colnames) - 1 else "")

	file = open("C:\\Users\\horat\\Desktop\\msc_2025\\R\\data\\busco_genes\\all_genes.tsv", "w")
	file.write(bin_genes_tsv)
	file.close()


if __name__ == "__main__":
	main()