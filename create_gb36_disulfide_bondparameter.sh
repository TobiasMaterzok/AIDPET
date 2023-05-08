#!/bin/bash

# 
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
# Date: 2020
# 
# This script modifies a GROMACS topology file to add
# "gb_36" to the end of bond lines that don't have it.
# It is designed to update disulfide bond parameters
# for the GROMOS 54a7 force field after cross-linking using pdb2gmx.
# In GROMACS <= 2021.1 pdb2gmx doesnt add bond coefficient to the
# disulfide bridge it generates into the topology file. Thus, we
# need to add them afterwards.
#

topology_file=$1

# Find the starting line number for the "bonds" section in the GROMACS topology file
bond_start=$(awk 'BEGIN{st=0}{if($2=="bonds") {print NR+2}}' $topology_file)

# Find the ending line number for the "bonds" section and the starting line number for the "pairs" section
bond_end=$(awk 'BEGIN{st=0}{if($2=="angles") {print NR-2}}' $topology_file)
pair_start=$(awk 'BEGIN{st=0}{if($2=="pairs") {print NR-2}}' $topology_file)

# Add "gb_36" to the end of bond lines that don't have it, and print the modified lines
awk -v "s=$bond_start" -v "e=$bond_end" -v "e2=$pair_start" '{
    if((NR>=s) && (NR<=e) && (NR<=e2)) {
        if(!index($4,"gb")) print $0,"gb_36";
        else print $0
    } else print $0
}' $topology_file
