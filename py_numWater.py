#!/usr/bin/env python

"""
Author: Tobias Materzok https://github.com/TobiasMaterzok/
Date: 2020

This script calculates the number of water molecules to add to a system to 
achieve a given weight fraction.

Usage: python py_numWater.py GRO_FILE WEIGHT_FRACTION

Input:
- GRO_FILE: the name of the Gromacs GRO file
- WEIGHT_FRACTION: the desired weight fraction of water to add

Output:
- The number of water molecules to add, printed to the console
"""

import MDAnalysis as md
import sys

# Read the command line arguments
gro = str(sys.argv[1])
wt = float(sys.argv[2])

# Create an MDAnalysis Universe object from the GRO file
u = md.Universe(gro)

# Select all atoms with a positive mass
sel = u.select_atoms('prop mass > 0.0')

# Calculate the total mass of the selected atoms
sel.accumulate('masses')
mass = sel.accumulate('masses')

# Calculate the weight of water to add based on the desired weight fraction
weightW = (wt * mass) / (1.0 - wt)

# Calculate the number of water molecules to add based on the weight of water and the mass of a single water molecule
num = int(round(weightW / 18))

# Print the number of water molecules to add to the console
print(num)
