import re
import os 
import json

def generate_tsv_from_busco_summaries(busco_summaries_folder):
	busco_summary_data = []

	for filename in os.listdir(busco_summaries_folder):
		json_file = open(busco_summaries_folder + "/" + filename)
		json_data = json.load(json_file)

		busco_data_entry = {}

		busco_data_entry["SampleName"] = re.findall("[R]+[0-9]+", json_data["parameters"]["out"])[0]
		busco_data_entry["Bin"] = re.findall("(?<=_000)[0-9]+", json_data["parameters"]["out"])[0]
		busco_data_entry["Specific"] = "1" if ("specific" in filename) else "0"

		busco_data_entry["LineageDataset"] = json_data["lineage_dataset"]["name"]
		busco_data_entry["LineageBUSCOs"] = json_data["lineage_dataset"]["number_of_buscos"]
		busco_data_entry["LineageSpecies"] = json_data["lineage_dataset"]["number_of_species"]

		busco_data_entry["Complete"] = json_data["results"]["Complete"]
		busco_data_entry["Single Copy"] = json_data["results"]["Single copy"]
		busco_data_entry["Multiple"] = json_data["results"]["Multi copy"]
		busco_data_entry["Fragmented"] = json_data["results"]["Fragmented"]
		busco_data_entry["Missing"] = json_data["results"]["Missing"]
		busco_data_entry["Number of Markers"] = json_data["results"]["n_markers"]
		busco_data_entry["Domain"] = json_data["results"]["domain"]
		busco_data_entry["Scaffolds"] = json_data["results"]["Number of scaffolds"]
		busco_data_entry["Contigs"] = json_data["results"]["Number of contigs"]
		busco_data_entry["Total Length"] = json_data["results"]["Total length"]
		busco_data_entry["Gaps"] = json_data["results"]["Percent gaps"]
		busco_data_entry["Scaffold N50"] = json_data["results"]["Scaffold N50"]
		busco_data_entry["Contig N50"] = json_data["results"]["Contigs N50"]

		busco_summary_data.append(busco_data_entry)


	entries_to_delete = []
	for i in range(len(busco_summary_data)):
		for j in range(len(busco_summary_data)):
			if busco_summary_data[i]["SampleName"] == busco_summary_data[j]["SampleName"] and busco_summary_data[i]["Bin"] == busco_summary_data[j]["Bin"]:
				print("Match: ", busco_summary_data[i]["Specific"], i, busco_summary_data[i]["Bin"], busco_summary_data[j]["Specific"], j,  busco_summary_data[j]["Bin"])
				if busco_summary_data[i]["Specific"] != busco_summary_data[j]["Specific"]:
					if busco_summary_data[i]["Specific"] == "0" and busco_summary_data[j]["Specific"] == "1":
						entries_to_delete.append(i)
					elif busco_summary_data[i]["Specific"] == "1" and busco_summary_data[j]["Specific"] == "0":
						entries_to_delete.append(j)
				else:
					if busco_summary_data[i] != busco_summary_data[j]:
						print("What!")
						print(busco_summary_data[i])
						print(busco_summary_data[j])
						print("========")

	print([busco_summary_data[idx]["Bin"] for idx in entries_to_delete])
	busco_summary_data = [busco_summary_data[idx] for idx in range(len(busco_summary_data)) if not(idx in entries_to_delete)]

	busco_summary_tsv = ""

	tsv_colnames = ["SampleName", "Bin", "Specific", "LineageDataset", "LineageBUSCOs", "LineageSpecies", "Complete", "Single Copy", "Multiple", "Fragmented", "Missing", "Number of Markers", "Domain", "Scaffolds", "Contigs", "Total Length", "Gaps", "Scaffold N50", "Contig N50"]

	for i in range(len(tsv_colnames)):
		busco_summary_tsv += tsv_colnames[i] + ("\t" if i < len(tsv_colnames) - 1 else "")

	for data_entry in busco_summary_data:
		busco_summary_tsv += "\n"
		for i in range(len(tsv_colnames)):
			busco_summary_tsv += str(data_entry[tsv_colnames[i]]) + ("\t" if i < len(tsv_colnames) - 1 else "")

	# I know, I know, it's hard-coded...!!
	file = open("C:\\Users\\horat\\Desktop\\msc_2025\\R\\data\\busco_tsv\\RAMAN_R018.tsv", "w")
	#file = open("../R/data/busco_tsv/RAMAN_R018.tsv", "w+")
	file.write(busco_summary_tsv)
	file.close()

	return busco_summary_tsv

if __name__ == "__main__":
	generate_tsv_from_busco_summaries("C:\\Users\\horat\\Desktop\\msc_2025\\R\\data\\busco_summaries\\RAMAN_R018")
	#generate_tsv_from_busco_summaries("../R/data/busco_summaries/RAMAN_R018")
