"""Contact Map Generator
Author-Geet Madhukar
Date 05/12/26"""



import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB

def get_ca_coordinate(pdb_file, chain_id='A'):
	parser=PDB.PDBParser(QUIET=True)
	structure=parser.get_structure('protein', pdb_file)

	coords=[]

	for residue in structure[0][chain_id]:
		if 'CA' in residue:
			coords.append(residue['CA'].get_vector().get_array())
	return np.array(coords)


def build_contact_map(coords):
	n=len(coords)
	m=np.zeros((n,n))
	for i in range(len(coords)):
		for j in range(len(coords)):
			dist=np.sqrt(np.sum((coords[i] - coords[j])**2))
			if dist <=8:
				m[i][j]=1

	return m

def main():
	import sys
	if len(sys.argv) != 3:
		print("Usage python contact_map pdb_file chain id")
		sys.exit(1)

	pdb_file=sys.argv[1]
	chain_id=sys.argv[2]

	coords=get_ca_coordinate(pdb_file, chain_id)

	contact_map=build_contact_map(coords)

	plt.figure(figsize=(8,8))
	plt.imshow(contact_map, cmap='Blues', origin='lower')
	plt.colorbar(label='Contact (1=yes, 0=no)')
	plt.xlabel('Residue Index')
	plt.ylabel('Residue Index')
	plt.title('Cα Contact Map')
	plt.tight_layout()
	plt.savefig('contact_map.png', dpi=150)
	print("Contact map saved as contact_map.png")

	
if __name__=="__main__":
		main()
