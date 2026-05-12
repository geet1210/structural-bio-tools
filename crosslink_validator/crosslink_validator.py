"""Crosslink Validator
Author-Geet Madhukar
Date 05/11/26"""


import pandas as pd
import numpy as np
from Bio import PDB

def load_crosslinks(csv_file):
	xl=pd.read_csv(csv_file)
	return xl

def validate_crosslinks(pdb_file, crosslinks, chain_id="A"):
	parser=PDB.PDBParser(QUIET=True)
	structure=parser.get_structure("protein", pdb_file)
	coords={}
	for residue in structure[0][chain_id]:
		if "CA" in residue:
			coords[residue.get_id()[1]] =residue["CA"].get_vector().get_array()

	for _, row in crosslinks.iterrows():
		res1 = int(row["residue1"])
		res2 = int(row["residue2"])
		max_dist = row["max_distance"]

		distance = np.sqrt(np.sum((coords[res1] - coords[res2])**2))

		if distance <= max_dist:
			print(f"Residues {res1}-{res2}: {distance:.2f}Å ✅ SATISFIED")
		else:
			print(f"Residues {res1}-{res2}: {distance:.2f}Å ❌ VIOLATED")

def main():
	import sys
	if len(sys.argv)!= 3:
		print("Usage script, pdb file, csv file")
		sys.exit(1)

	pdb=sys.argv[1]
	csv_file=sys.argv[2]
	

	crosslinks = load_crosslinks(csv_file)
	validate_crosslinks(pdb, crosslinks)

if __name__=="__main__":
	main()

	


