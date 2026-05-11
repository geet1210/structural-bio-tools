"""RMSD Calculator
Author-Geet Madhukar
Date 05/11/26"""

import numpy as np
from Bio import PDB

def get_ca_coordinate(pdb_file, chain_id='A'):
	parser=PDB.PDBParser(QUIET=True)
	structure=parser.get_structure("protein",pdb_file)
	coords=[]
	for residue in structure[0][chain_id]:
		if "CA" in residue:
			coords.append(residue["CA"].get_vector().get_array())
	return np.array(coords)

def cal_rmsd(coords1, coords2):
	if len(coords1)!= len(coords2):
		raise ValueError("Structure diff lengths")
	diff =coords1-coords2

	return np.sqrt(np.mean(np.sum(diff**2,axis=1)))

def main():
	import sys
	if len(sys.argv) != 4:
		print("Usage python rmsd.py structure1.pdb structure2.pdb chain_id")
		sys.exit(1)

	pdb1=sys.argv[1]
	pdb2=sys.argv[2]
	chain=sys.argv[3]

	coords1=get_ca_coordinate(pdb1, chain)
	coords2=get_ca_coordinate(pdb2, chain)

	rmsd=cal_rmsd(coords1, coords2)
	print(f"RMSD between {pdb1} and {pdb2} = {rmsd:.3f} Å")

if __name__ == "__main__":
    main()