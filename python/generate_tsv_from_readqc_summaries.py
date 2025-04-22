import re
import os 

def generate_tsv_from_readqc_summaries(run_accessions, run_directory=""):
	#print(run_accessions)

	readqc_summary_data = {
#		"SRR29980928" : {
#			"pre" : {
#				"number_of_reads" : 9999,
#				"mean_read_length" : 9999,
#				"read_n50" : 9999,
#				"mean_read_quality" : 9999,
#				"total_bases_mbp" : 9999,
#			}
#		}
		 
	}

	for run_acc in run_accessions:
		summary_file = open((run_directory) + run_acc + "/summary/" + run_acc + ".RUN01.readqc_report.html")
		summary_file_lines = summary_file.readlines()

		run_acc_data = {}

		for i in range(len(summary_file_lines)):
			if "Number of reads" in summary_file_lines[i]:
				run_acc_data["number_of_reads"] = summary_file_lines[i + 1].strip().replace("<td>", "").replace("</td>", "")
			if "Mean read length" in summary_file_lines[i]:
				run_acc_data["mean_read_length"] = summary_file_lines[i + 1].strip().replace("<td>", "").replace("</td>", "")
			if "<b>Read N" in summary_file_lines[i]:
				run_acc_data["read_n50"] = summary_file_lines[i + 1].strip().replace("<td>", "").replace("</td>", "")
			if "Mean read quality" in summary_file_lines[i]:
				run_acc_data["mean_read_quality"] = summary_file_lines[i + 1].strip().replace("<td>", "").replace("</td>", "")
			if "Total bases" in summary_file_lines[i]:
				run_acc_data["total_bases_mbp"] = summary_file_lines[i + 1].strip().replace("<td>", "").replace("</td>", "")

		contig_coverage_loc = (run_directory) + run_acc + "/assembly/contig_qc/coverage/" + run_acc + ".coverage.txt"

		if os.path.exists(contig_coverage_loc):
			contig_coverage_file = open(contig_coverage_loc)
			contig_coverage_lines = contig_coverage_file.readlines()
			run_acc_data["contig_count"] = str(len(contig_coverage_lines) - 1)

		else:
			run_acc_data["contig_count"] = "0"


		assembly_stats_loc = (run_directory) + run_acc + "/assembly/bin_QC/assembly_stats/" + run_acc + ".RUN01.assembly_stats.csv"

		if os.path.exists(assembly_stats_loc):
			assembly_stats_file = open(assembly_stats_loc)
			assembly_stats_lines = assembly_stats_file.readlines()
			number_of_bins = len(assembly_stats_lines) - 1

			run_acc_data["bin_count"] = str(number_of_bins)
			
			t = 0
			for i in range(len(assembly_stats_lines)):
				if i != 0:
					split_line = assembly_stats_lines[i].split(",")
					t += int(split_line[7].strip())

			run_acc_data["average_bin_n50"] = str(float(t) / float(number_of_bins))

		else:
			run_acc_data["bin_count"] = "0"
			run_acc_data["average_bin_n50"] = "0"

		readqc_summary_data[run_acc] = run_acc_data

	readqc_summary_tsv = "run_accession	"
	readqc_summary_colnames = ["number_of_reads", "mean_read_length", "read_n50", "mean_read_quality", "total_bases_mbp", "contig_count", "bin_count", "average_bin_n50"]

	for i in range(len(readqc_summary_colnames)):
		readqc_summary_tsv += readqc_summary_colnames[i] + ("\t" if (i + 1) < len(readqc_summary_colnames) else "")

	for run_acc in readqc_summary_data.keys():
		readqc_summary_tsv += "\n" + run_acc + "\t"
		for i in range(len(readqc_summary_colnames)):
			readqc_summary_tsv += readqc_summary_data[run_acc][readqc_summary_colnames[i]] + ("\t" if (i + 1) < len(readqc_summary_colnames) else "")

	return readqc_summary_tsv
