#!/usr/bin/python

"""
	The script cross checks that species names in old alignment files and old 
	Rdata_file correspond.  
"""

a = []
b = []

with open("/mnt/c/WB/tRNA/_old/R_data.txt", "rt") as infile:
	for line in infile.readlines():
		line = line.split()
		a.append(line[0])

with open("/mnt/c/WB/tRNA/_old/020_tetrapodes_tRNA/3_FINAL_ALIGNMENT/A.tRNA", "rt") as infile:
	for line in infile.readlines():
		line = line.split()
		b.append(line[0])

#check if overlap generator approach 
print(set(a) & set(b))
print(len(a), len(b))
print(len(set(a) & set(b)))