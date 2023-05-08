"""
Author: Tobias Materzok https://github.com/TobiasMaterzok/
Date: 2021

This script reads a protein GROMACS file (.gro) and counts the number of occurrences
of a specified atom type in the system. The script takes two command-line arguments:
the path to the protein GROMACS file (.gro) and the atom type (e.g., "SG" for the 
sulfur atom which is specific for the cysteine residue in the GROMOS 54A7 FF)
to count in the system. It then reads the coordinates of the specified atom type and
prints the count of those atoms in the first frame of the trajectory.

Input:
- traj_file: Path to the protein GROMACS file (.gro)
- atom_type: Atom type to count in the system (e.g., "SG")

Output:
- Prints the count of specified atom type in the first frame of the trajectory
"""

import numpy as np
import sys

# Input arguments
traj_file = sys.argv[1]
atom_type = str(sys.argv[2])

# Initialize variables
coordinates = []
box_dimensions = []
frame_count = 0

# Read coordinates
with open(traj_file, 'r') as file:
    for line in file:
        coordinates.append([])
        box_dimensions.append([])
        natoms = int(file.readline().strip())
        
        for i in range(natoms):
            row = file.readline()
            current_atom_type = row[10:15].strip()
            
            if current_atom_type == atom_type:
                coords = np.array([row[20:28].split(), row[28:35].split(), row[36:44].split()]).astype(float)
                coordinates[frame_count].append(coords)
        
        box_dimensions[frame_count].append(np.array(file.readline().split()).astype(float))
        frame_count += 1

coordinates = np.array(coordinates)
box_dimensions = np.array(box_dimensions)

selected_atoms = len(coordinates[0])
print(selected_atoms)
