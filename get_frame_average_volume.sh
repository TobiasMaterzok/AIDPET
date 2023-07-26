#!/bin/bash

#
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
# Date: 2021/2022
#
# This script processes a GROMACS trajectory file (XTC) and energy file (EDR) to extract the frame with 
# the volume value closest to the average volume. The average volume is calculated from the energy file 
# after a specified start point. Only frames at intervals specified by the energy output interval 
# are considered. The selected frame is then extracted and saved as a GRO file.
# 
# Pressure fluctuations in molecular dynamics simulations can be high due to the dynamic nature of 
# interatomic forces, which cause the instantaneous pressure values to vary significantly during 
# the simulation. This is because the virial pressure calculation in molecular dynamics depends on the
# instantaneous positions and forces of all particles in the system, which are constantly changing as the 
# simulation progresses. The root mean square deviation of the pressure for systems like dry gecko keratin can 
# reach around 300 to 400 bar.
# 
# Focusing on instantaneous pressure values might lead to selecting a non-representative configuration. 
# Instead, this script selects a configuration where the volume is closest to the volume average of 
# the last ~20% of the NPT equilibration. Volume, as a thermodynamic property, is less affected by 
# pressure fluctuations and provides a more stable reference for selecting a representative configuration.
# 
# By choosing a frame with a volume close to the average volume, we ensure that the configuration has 
# experienced a sufficiently long equilibration time and is representative of the equilibrium state.
# 
# Input:
# 
#     deffnm - The basename of the GROMACS run files (e.g., 'npt' if files are named 'npt.edr', 'npt.xtc', etc.)
#     trjoutinterval_by_energyoutinterval - The ratio of trajectory output interval to energy output interval
# 
# Output:
# 
#     confout_avgvol.gro - A configuration file with the closest volume to the average 
#     volume of the last ~20% of the NPT equilibration
#

tools=~/AIDPET
. $tools/functions.sh
load_gromacs_2018_4

# Set input arguments and load necessary modules
deffnm=$1
trjoutinterval_by_energyoutinterval=$2 # nstxout-compressed / nstcalcenergy = 10
tbe=$trjoutinterval_by_energyoutinterval

# Extract volume information from the .edr file
volume_gmxname=$(~/tools/get_number_for_gmxenergy.sh "$deffnm".edr Volume)
echo "$volume_gmxname " | gmx energy -f "$deffnm".edr -o gmx_volume.xvg

# Determine the starting and end point for the analysis
end_line=$(cat gmx_volume.xvg | awk '{if($1=="@") print NR+1}' | tail -1)
tail -n +"$end_line" gmx_volume.xvg > tmptmptmp && mv tmptmptmp gmx_volume.xvg
start_point=$(tail -1 gmx_volume.xvg | awk '{print $1/1.2}')
echo $start_point

# Calculate the average volume
echo "$volume_gmxname " | gmx energy -f "$deffnm".edr -o tmptmp.xvg -b "$start_point" 1> voltmp.xvg
avgvol=$(tail -1 voltmp.xvg | awk '{print $2}')
cat voltmp.xvg
rm voltmp.xvg tmptmp.xvg

# Find the frame with the volume closest to the average volume, considering frames after the start point and at energy_out_interval
# Print the frame's time value, volume value, and the squared difference between the volume and the average volume, then sort by the squared difference
# Assign the time value of the frame with the smallest squared difference to the 'selected_frame' variable
selected_frame=$(awk -v "s=$start_point" -v "step=$tbe" -v "volavg=$avgvol" '{if($1>s){if((volavg-$2)^2 < 1.0){if($1%step == 0) print $1,$2,(volavg-$2)^2  } }}' gmx_volume.xvg | sort -k3 | head -1 | awk '{print $1}')
rm gmx_volume.xvg

# Extract the frame with the closest volume to the average and save it as a .gro file
echo "0 " | gmx trjconv -f "$deffnm".xtc -s "$deffnm".tpr -o confout_avgvol.gro -b "$selected_frame" -e "$selected_frame"
