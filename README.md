# Structural Biology Tools
A collection of modular, reusable computational structural biology 
algorithms written in Python from scratch.

## Tools

### 1. RMSD Calculator
Calculate RMSD between two protein structures using Cα atoms.
- **Input:** Two PDB files, chain ID
- **Usage:** `python rmsd_cal.py structure1.pdb structure2.pdb chain_id`

### 2. Cross-link Validator
Validate XL-MS crosslinks against a PDB structure.
- **Input:** PDB file, CSV of crosslink pairs
- **Usage:** `python crosslink_validator.py structure.pdb crosslinks.csv`

### 3. Contact Map Generator
Generate and visualize Cα contact maps from PDB structures.
- **Input:** PDB file, chain ID
- **Usage:** `python contact_map.py structure.pdb chain_id`

## Requirements
```bash
pip install biopython numpy pandas matplotlib
```

## Author
[geet1210](https://github.com/geet1210) — Computational Structural Biologist