# 
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
# Date: 2020
#
# This script generates a GROMACS .gro and .top file for a given input .pdb file for the GROMOS54A7
# force field, disulfide bonds, termini, and other interactive options. It automates the 
# process of building the molecular system for GROMACS simulations, allowing users to easily generate 
# the necessary files for their specific system.
# 
# Input:
# 
#     pdb_input: The input protein structure in .pdb format (hardcoded as "all.pdb")
# 
# Output:
# 
#     gro_output: The output GROMACS .gro file containing the system's coordinates (provided as the first argument)
#     top_output: The output GROMACS .top file containing the system's topology (provided as the second argument)
# 
# The script performs the following steps:
# 
#     Runs pdb2gmx to select amino acid types for non-standard residues.
#     Runs a second pdb2gmx command to select all cysteines for disulfide bond formation.
#     Runs a third pdb2gmx command to select the termini of the protein.
#     Combines the selections from steps 1-3 and runs the final pdb2gmx command with interactive mode 
# 	and ignoring hydrogen atoms.
# 
# Usage:
# script.sh gro_output top_output
# 

tools=~/tools_ua_gecko
source ~/tools_ua_gecko/functions.sh
load_gromacs_2018_4

gro_output="$1"
top_output="$2"
pdb_input="all.pdb"

force_field="gromos54a7"

# Run pdb2gmx to select amino acid types for non-standard residues
yes "0" | gmx pdb2gmx -ff "$force_field" -f "$pdb_input" -o "$gro_output" -water spce -p "$top_output" -merge all -chainsep ter -lys -arg -asp -glu -gln > sel_CA.dat
sel_CA=$(cat sel_CA.dat | grep "charge 0" | awk '{printf("%s\\n",substr($1,1,1))}')
rm "$gro_output".gro sel_CA.dat "$top_output".top posre.tip

# Run the second pdb2gmx command to select all cysteines
yes "n" | gmx pdb2gmx -ff "$force_field" -f "$pdb_input" -o "$gro_output" -water spce -p "$top_output" -merge all -chainsep ter -ss 2> sel_SS.dat
sel_SS=$(cat sel_SS.dat | grep -o -i Link | awk '{printf("n\\n")}')
rm "$gro_output".gro sel_SS.dat "$top_output".top posre.tip

# Run the third pdb2gmx command to select termini
yes "0" | gmx pdb2gmx -ff "$force_field" -f "$pdb_input" -o "$gro_output" -water spce -p "$top_output" -merge all -chainsep ter -ter 1> sel_termini.dat
sel_termini=$(cat sel_termini.dat | grep Select | awk '{printf("0\\n")}')
rm "$gro_output".gro sel_termini.dat "$top_output".top posre.tip

# Combine selections and run the final pdb2gmx command with interactive mode and ignoring hydrogen atoms
str="$sel_CA""$sel_SS""$sel_termini"

printf $str | gmx pdb2gmx -ff "$force_field" -f "$pdb_input" -o "$gro_output" -water spce -p "$top_output" -merge all -chainsep ter -inter -ignh
