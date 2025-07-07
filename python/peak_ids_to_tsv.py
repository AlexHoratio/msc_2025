import os
import re

###################################
#
# This script will read in a series of Peak ID summaries from the "Peak IDs" folder (must be placed next to this script)
# 	and will create a folder "Peak ID TSVs" full of formatted TSV files containined peak RTs, Area %, three reported identifications, and quality for each
#
###################################

def main():
	re.DOTALL = True
	for filename in os.listdir("Peak IDs"):
		peak_ids_tsv = "Peak\tRT\tArea (%)\tID #1\tRef #1\tCAS #1\tQual #1\tID #2\tRef #2\tCAS #2\tQual #2\tID #3\tRef #3\tCAS #3\tQual #3"

		peak_ids_txt = open("Peak IDs/" + filename, "r")
		peak_id_line = ""

		making_compound_name = False
		compound = ""
		compound_name = ""
		for line in peak_ids_txt.readlines()[18::]:
			if re.search("\\A[\\s]{0,2}[0-9]", line): #is "header line" for a peak
				if peak_id_line != "": # if the previous line was something, it is now done, write it to the tsv
					if making_compound_name: #last compound name is finished
						compound = compound_name + compound
						
						peak_id_line += "\t" + compound
						compound_name = ""
						making_compound_name = False

					peak_ids_tsv += "\n" + peak_id_line
					peak_id_line = ""

				header = [x for x in line.replace("\n", "").split(" ") if x != ""]
				peak_id_line += header[0] + "\t" + header[1] + "\t" + header[2]

			if re.search("\\A[\\s]{17}", line): # is any line containing compound title
				compound_line = [x for x in line.replace("\n", "").split(" ") if x != ""]

				if re.search(".{72}[0-9]{1,2}[\\n]\\Z", line): # is the first line for the compound name
					if making_compound_name: #last compound name is finished
						compound = compound_name + compound
						
						peak_id_line += "\t" + compound
						compound_name = ""
						making_compound_name = False

					making_compound_name = True
					
					compound = "\t" + compound_line[-3] + "\t" + compound_line[-2] + "\t" + compound_line[-1]

				if making_compound_name:
					compound_name = compound_name + " ".join([x for x in compound_line if re.search("[a-zA-Z]+", x)])
					
		if not(os.path.exists("Peak ID TSVs")):
			os.makedirs("Peak ID TSVs")

		output = open("Peak ID TSVs/" + filename, "w")
		output.write(peak_ids_tsv)
		output.close()

if __name__ == "__main__":
	main()