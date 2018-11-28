#!/usr/bin/python

"""The script takes old alignment files from 020_tetrapodes_tRNA, for each
	tRNA in each animals it records into a dict of dicts the number of paired bases in each tRNA and number of paired
	bases that are either G or C.

	After collection the info is appended to the R_data.txt file"""

import os
import re
import numpy # for division by zero later on

path = '/mnt/c/WB/tRNA/_old/020_tetrapodes_tRNA/3_FINAL_ALIGNMENT/'

# List of files with .tRNA extension in path
ls = filter(lambda x: x.endswith(('.tRNA')), os.listdir(path))
#print(ls)

# Amino-acid abbreviations dict
aa_abbr = {
			'A':'Ala',
			'C':'Cys', 
			'D':'Asp', 
			'E':'Glu', 
			'F':'Phe', 
			'G':'Gly', 
			'H':'His', 
			'I':'Ile', 
			'K':'Lys', 
			'L1':'LeuUUR',
			'L2':'LeuCUN', 
			'M':'Met', 
			'N':'Asn', 
			'P':'Pro', 
			'Q':'Gln', 
			'R':'Arg', 
			'S1':'SerUCN', 
			'S2':'SerAGY', 
			'T':'Thr', 
			'V':'Val', 
			'W':'Trp', 
			'Y':'Tyr'}

# THERE MUST BE A MORE ELEGANT WAY TO DO THIS! BUT FOR NOW
# This bit of code gets a list of species, creates the whacky dict
stats = {}

with open("/mnt/c/WB/tRNA/_old/020_tetrapodes_tRNA/3_FINAL_ALIGNMENT/A.tRNA", "rt") as infile:
	for line in infile.readlines():
		line = line.split()
		stats[line[0]] = {}
		for i in aa_abbr:
			stats[line[0]][aa_abbr[i]] = {}		

# Collecting the data from old alignment files and storing them in stats
for x in ls:
	aa = aa_abbr[os.path.splitext(x)[0]] # Filename wo an extension
	#print(iks)
	with open("%s%s" % (path, x), "rt") as infile:
		for line in infile.readlines():
			line = line.split()
			string = line[-1]
			#print(line[0])
			#if WC == 0: # F*ing Marsupials!
			
			lines = [m.start() for m in re.finditer('\|', string)]
			trna = string[lines[0]+1:lines[-1]]
			trna = re.sub('-', '', trna)
			trna = re.sub('\|', '', trna)
			WC = (sum(1 for letter in trna if letter.isupper()))
			GC = (sum(1 for c in trna if c == "G" or c == "C"))
			stats[line[0]][aa]["WC"] = WC
			stats[line[0]][aa]["GC"] = GC
			stats[line[0]][aa]["len"] = len(trna)
			print(trna)
			print(WC,GC)
			#print(WC, GC)
# Data is gathered. Now we open old R_data file, create new one and dump data
outfile = open("/mnt/c/WB/tRNA/2_processed_files/5_from_old/R_data_1_WCGC.txt","wt")
with open("/mnt/c/WB/tRNA/_old/R_data.txt", "rt") as infile:
	first_line = True
	for line in infile.readlines():
		line = line.split()
		if "sister_pairs" in line:
			break
		if first_line == True:
			for i in aa_abbr:
				line.append("%s_len" % aa_abbr[i])
				line.append("%s_WC" % aa_abbr[i])
				line.append("%s_GC" % aa_abbr[i])
				first_line = False
		else:
			for i in aa_abbr:
				line.append(str(stats[line[0]][aa]["len"]))
				line.append(str(stats[line[0]][aa_abbr[i]]["WC"]))
				line.append(str(stats[line[0]][aa_abbr[i]]["GC"]))
		outfile.write('\t'.join(line) + '\n')		
outfile.close()