#!/bin/bash

#
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
# Date: 2020
# 
# This script performs a step-by-step process to generate a GROMACS .gro and .top file for a given 
# input .pdb file, considering the GROMOS54A7 force field, disulfide bonds, termini, and other interactive options. 
# It automates the process of building the molecular system for GROMACS simulations, allowing users to 
# easily generate the necessary files for their specific system.
#
# Input:
#
#    input_pdb_file: The input protein structure in .pdb format
#
# Output:
#
#    output_gro_file: The output GROMACS .gro file containing the system's coordinates
#    output_top_file: The output GROMACS .top file containing the system's topology
#
# Usage:
# script.sh output_gro_file output_top_file input_pdb_file
#

gro_output="$1"
top_output="$2"
pdb_input="$3"

tools=~/tools_ua_gecko
source ~/tools_ua_gecko/functions.sh
load_gromacs_2018.4

force_field="gromos54a7"

# Run pdb2gmx to select amino acid types for non-standard residues
yes "0" | gmx pdb2gmx -ff "$force_field" -f "$pdb_input" -o "$gro_output" -water spce -p "$top_output" -merge all -chainsep ter -lys -arg -asp -glu -gln > sel_CA.dat
sel_CA=$(cat sel_CA.dat | grep "charge 0" | awk '{printf("%s\\n",substr($1,1,1))}')
rm "$gro_output".gro sel_CA.dat "$top_output".top posre.tip

# Run pdb2gmx to identify disulfide bonds
gmx pdb2gmx -ff "$force_field" -f "$pdb_input" -o "$gro_output" -merge all -water spc -chainsep ter -ignh 2> sel_SS.dat
cat sel_SS.dat | grep Linking > sel_SS_tmp.dat
# Identify unique disulfide bonds involving cysteine residues
# Output the number of unique disulfide bonds found in the input file
sel_SS=$(cat sel_SS_tmp.dat | awk '{a=substr($3,4); b=substr(substr($6,1,length($6)-3),4); bA=1; bB=1; for (i=0;i<NR;i++){if(a==X[i] || a==Y[i] || b==X[i] || b==Y[i]){bA=0}} if(bA==1){X[NR]=a; Y[NR]=b}; Accept="y"; for (i=0;i<NR;i++){if(X[NR]==X[i] || Y[NR]==X[i] || X[NR]==Y[i] || Y[NR]==Y[i]) {Accept="n"}}; printf("%s\\n",Accept)}')
rm "$gro_output".gro sel_SS.dat "$top_output".top posre.tip

# Run pdb2gmx to select termini types
yes "0" | gmx pdb2gmx -ff "$force_field" -f "$pdb_input" -o "$gro_output" -water spce -p "$top_output" -merge all -chainsep ter -ter 1> sel_termini.dat
sel_termini=$(cat sel_termini.dat | grep Select | awk '{printf("0\\n")}')
rm "$gro_output".gro sel_termini.dat "$top_output".top posre.tip

# Combine the selections and rerun pdb2gmx with interactive options
str="$sel_CA""$sel_SS""$sel_termini"
printf $str | gmx pdb2gmx -ff "$force_field" -f "$pdb_input" -o "$gro_output" -water spce -p "$top_output" -merge all -chainsep ter -inter -ignh -ss

# Adjust topology for specific force field (gromos54a7)
$tools/create_gb36_disulfide_bondparameter.sh "$top_output" > tmptop && mv tmptop "$top_output"
