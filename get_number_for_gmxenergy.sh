#!/bin/bash

#
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
# Date: 2019
# 
# This script finds and prints the index number of a
# specified energy term in a GROMACS energy file.
# 

# Assign command line arguments to variables
energy_file_path=$1
energy_term=$2

source ~/tools_ua_gecko/functions.sh
load_gromacs_2018.4

# Run the gmx energy command and create temporary files for output and error streams
echo " " | gmx energy -f "$energy_file_path" -o energy_values_tmp 2> energy_terms_list_tmp

# Extract the term from the error stream (energy_terms_list_tmp) using grep and awk
cat energy_terms_list_tmp | grep -o -P ".{0,5}$energy_term " | awk '{print $1}' | head -1

# Clean up the temporary files created during the process
rm energy_terms_list_tmp
