#!/bin/bash

#
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
# Date: 2020
#
# Description:
# This script calculates the number of disulfide bonds required to achieve a target percentage of disulfide bond density,
# then iteratively increases the cutoff distance for disulfide bonds and runs GROMACS pdb2gmx to find disulfide bonds
# within that cutoff distance. The process continues until the target number of disulfide bonds is reached or the maximum
# number of iterations (60) is exceeded.
#
# Inputs:
#   1. deffnm: Prefix for the GROMACS files (e.g., "my_simulation").
#   2. num_atoms_in_chain: Number of atoms in each protein chain.
#   3. percentage_ssdensity: Target percentage of cysteins being involved in a disulfide bond (0-1, e.g., 0.33 for 33%).
#
# Outputs:
#   - GA9rnd_gromos_ss.gro and GA9rnd_gromos_ss.top: GROMACS structure and topology files with the selected disulfide bonds.
#   - ERROR file with error messages if the process fails (e.g., not finding enough disulfide bonds or missing trajectory file).
#

deffnm="$1"
num_atoms_in_chain=$2
percentage_ssdensity=$3
lastatomAfterGMXTrjconv=O2
naic=$num_atoms_in_chain
tools=~/AIDPET
source $tools/functions.sh
load_gromacs_2018_4


force_field=gromos54a7

# Create a new directory for ssbond and move to the parent directory
mkdir ssbond
cd ssbond 
cd ..
# Extract the last frame from the trajectory without hydrogen atoms
echo "0 " | gmx trjconv -f "$deffnm".xtc -s "$deffnm".tpr -o tmp.gro -dump 10000000000 2> lastframe
lastframe=$(cat lastframe | grep "last frame" | tail -1 | awk '{print substr($8,3,1000)}') # find the last frame in the trajectory
echo "2 " | gmx trjconv -f "$deffnm".xtc -s "$deffnm".tpr -o confout_pureP.pdb -dump $lastframe -pbc mol  # output only Protein-H
echo "2 " | gmx trjconv -f "$deffnm".xtc -s "$deffnm".tpr -o tmp.gro -dump $lastframe  # it needs to be without H atoms.. so output only Protein-H

# Move to the ssbond directory and copy the required files
cd ssbond
cp ../confout_pureP.pdb .
mv ../tmp.gro .

# Separate individual protein chains in the confout_pureP.pdb file
cat confout_pureP.pdb | awk -v "naic=$naic" -v "la=$lastatomAfterGMXTrjconv" '{if(($5%naic == 0) && ($3==la)){printf("%s\nTER\n",$0)} else print $0}' > sep.pdb

# Count the number of sulfur atoms (SG) in the system
N_SG=$(python $tools/py_get_number_of_atomtype.py tmp.gro "SG")

# Calculate the number of disulfide bonds needed to achieve the target percentage
N_bonds_needed=$(python -c "print($N_SG*$percentage_ssdensity)")
N_bonds_needed=$(echo "($N_bonds_needed/1)" | bc)

# Initialize variables
l=0
N_SS=0
count=0

# Check if the trajectory file exists
if [ -f ../"$deffnm".xtc ]
then
 # Loop until enough disulfide bonds are found or the limit of 60 iterations is reached (which is ~ 3nm radius)
while (( $N_SS < $N_bonds_needed ));
do
	((count++))
	if (( $count > 60 ))
	then
		echo "FAILED cant find enough ssbonds" >> ../ERROR
		exit 1
	fi

	# Update the cutoff value for disulfide bonds
	l=$(echo "$l + 0.05" | bc -l)
	cutoff=$(echo "0.2+$l" | bc -l)

	# Prepare the specbond.dat file
	mindist=0.1
	n=$(echo "($cutoff-$mindist)/0.01" | bc -l)
	N=$(echo "($n/1) + 1" | bc)
	echo "$N" > specbond.dat
	for i in `seq $mindist 0.01 $cutoff`
	do
		echo "CYS     SG      1       CYS     SG      1       $i     CYS2    CYS2" >> specbond.dat
	done
	
	# Run GROMACS pdb2gmx with the updated specbond.dat file
	gmx pdb2gmx -ff "$force_field" -f sep.pdb -o tmp.gro -merge all -water spc -chainsep ter -ignh 2> sel_SS.dat
	cat sel_SS.dat | grep Linking > sel_SS_tmp.dat

	# Count the number of unique disulfide bonds by comparing residue indices and ensuring no shared residues between bonds
	# Process sel_SS_tmp.dat using awk to extract residue indices of disulfide bonds, store unique bonds in X and Y arrays, and count the number of accepted (unique) bonds, ensuring that no residues are shared between different bonds
	N_SS=$(cat sel_SS_tmp.dat | awk '{a=substr($3,4); b=substr(substr($6,1,length($6)-3),4); bA=1; bB=1; for (i=0;i<NR;i++){if(a==X[i] || a==Y[i] || b==X[i] || b==Y[i]){bA=0}} if(bA==1){X[NR]=a; Y[NR]=b}; Accept="y"; for (i=0;i<NR;i++){if(X[NR]==X[i] || Y[NR]==X[i] || X[NR]==Y[i] || Y[NR]==Y[i]) {Accept="n"}}; printf("%s\n",Accept)}'| grep y | wc -l)
	
	# Clean up temporary files
	rm tmp.gro topol.top posre.tip

	# Print the current cutoff value and the number of disulfide bonds found
	echo "$cutoff $N_SS"

	# Wait for 1 second before the next iteration
	sleep 1
done

# Finalize the structure with the selected disulfide bonds
$tools/select_interactive_pdb2gmx_for_ssbonds_lightH.sh GA9rnd_gromos_ss.gro GA9rnd_gromos_ss.top sep.pdb
else
	echo "FAILED.. in create_amorph_ssbonds.sh missing "$deffnm".xtc trajectory" >> ../ERROR
fi
