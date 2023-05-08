#!/bin/bash

#
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
#
# A small bash script that combines useful functions for daily usage on a HPC Cluster
#

tools=~/tools_ua_gecko

function load_gromacs_2018.4(){
    module purge && module load intel/env/2018 fftw/intel/single/sse/3.3.8 gromacs/nompi/cpu/intel/single/2018.4
}

function load_gromacs_2021.1(){
    module purge && module load gcc/10.3.0 fftw/gcc/single/sse/3.3.9 gromacs/nompi/cpu/gcc/single/2021.1
}


# This function creates a template file for the Batch partition with specified parameters.
# Usage: batch_template <filename> <num_cpus> <string1> <string2>
# Arguments:
#   <filename>    : Output filename for the batch template.
#   <num_cpus>    : Number of CPUs to be used.
#   <string1>     : Additional user-specified string.
#   <string2>     : Another additional user-specified string.
function batch_template(){
        echo "#NCPU $2" >> $1
        echo "#TIME 4-00:00:00" >> $1
        echo "#$3" >> $1
        echo "#$4" >> $1
        echo "#" >> $1
}

# This function creates a template file for the Haswell partition with specified parameters.
# Usage: haswell_template <filename> <num_cpus> <string1> <string2>
# Arguments:
#   <filename>    : Output filename for the batch template.
#   <num_cpus>    : Number of CPUs to be used.
#   <string1>     : Additional user-specified string.
#   <string2>     : Another additional user-specified string.
function haswell_template(){
        echo "#NCPU $2" >> $1
        echo "#TIME 4-00:00:00" >> $1
        echo "#HASWELL 1" >> $1
        echo "#$3" >> $1
        echo "#$4" >> $1
}

# This function replaces a specified pattern in an input file with a new value.
# Usage: change_mdp <pattern> <replacement> <input_file>
# Arguments:
#   <pattern>     : The pattern to search for in the input file.
#   <replacement> : The replacement value for the pattern.
#   <input_file>  : The input file to modify.
function change_mdp(){
awk -v "pttrn=$1" -v "replace=$2" '{if($1 == pttrn){printf("%s             = %s\n",pttrn,replace)}else print $0}' $3 > tmp && mv tmp $3
}

# This function calculates the number of steps required for a GROMACS simulation.
# Usage: ns_to_steps <simulation_time> <time_step>
# Arguments:
#   <simulation_time> : Desired simulation time in nanoseconds.
#   <time_step>       : Time step size for the simulation in fs (femtoseconds).
function ns_to_steps(){
python -c "print(int($1*1000/$2))"
}


function awkprint1(){
awk '{printf("%s ",$1)}'
}

function awkprint2(){
awk '{printf("%s ",$2)}'
}


# This function finds the closest value in a specified column to a given target value. 
# Then, it prints the corresponding value from another specified column.
# Arguments:
# $1 (inputFile): The file to work on.
# $2 (targetColumn): The column to find the closest value in (columns are counted starting with 1).
# $3 (targetValue): The value to find the closest value to.
# $4 (outputColumn): The column to print the corresponding value from once the closest value is found.
#
# Example:
# test.xvg
# 11.100000 36.829657 6.101642
# 11.400000 36.402069 5.731998
# 11.700000 35.953025 5.372652
# closest_value_col_val test.xvg 2 36.399999 2
#   outputs: 36.402069
# closest_value_col_val test.xvg 2 36.399999 3
#   outputs: 5.731998
function closest_value_col_val(){
exa=$(cat $1 | awk -v c=$2 -v t=$3 'NR==1{d=$c-t;d=d<0?-d:d;v=$c;next}{m=$c-t;m=m<0?-m:m}m<d{d=m;v=$c}END{print v}')
cat $1 | awk -v c=$2 -v c2=$4 -v e=$exa '{if($c==e) print $c2}'
}


# This function modifies the job submission file for an HPC cluster.
# Usage: change_sbatch <input_file>
# Arguments:
#   <input_file> : The input file containing the parameters for the GROMACS simulation.
#
# Description:
#   The function reads the input file and extracts the required parameters for the GROMACS simulation,
#   such as the number of CPUs, runtime, GPU usage, and other settings. Based on these parameters,
#   it then updates the job submission file with the appropriate SLURM directives, modules, and
#   GROMACS commands. This makes it easier to submit GROMACS simulations to an HPC cluster.
#
# Input file format:
#   The input file should contain the following parameters:
#   #NCPU: Number of CPUs to use for the simulation
#   #TIME: The maximum runtime for the job in the format days-hours:minutes:seconds (e.g., 3-00:00:00)
#   #GPU: Set to 1 if using a GPU for the simulation, otherwise 0
#   #HASWELL: Set to 1 if using Intel Haswell architecture, otherwise 0
#   #2021: Set to 1 if using GROMACS 2021 version, otherwise 0
#   #GROUPCUTOFF: Set to 1 if using a group cutoff scheme, otherwise 0
#
# Example:
#   To modify a job submission file named "input.txt":
#   change_sbatch input.txt
function change_sbatch() {
        # Find the line number with the last occurrence of "##########################"
        lines=$(awk '{if($1=="##########################") print NR}' $1 | tail -1)

        # If the line number is found, remove all the lines above it (including the line itself) from the input file
        if (( $lines ))
        then
                lines=$(( lines+2 ))
                tail -n +"$lines" $1 > SBTMP && mv SBTMP "$1"
        fi

        # Extract relevant parameters from the input file
        NCPU=$(head -5 "$1" | awk '{if($1=="#NCPU") print $2}')
        if (( NCPU == 1 ))
        then
                part=serial
        else
                part=batch
        fi
        TIME=$(head -5 "$1" | awk '{if($1=="#TIME") print $2}')
        GPU=$(head -5 "$1" | awk '{if($1=="#GPU") print $2}')
        HASWELL=$(head -5 "$1" | awk '{if($1=="#HASWELL") print $2}')
        TWOTHOUSANDTWENTYONE=$(head -5 "$1" | awk '{if($1=="#2021") print $2}')

        # Determine the partition type based on the input parameters
        if (( GPU == 1 ))
        then
                part=batchgpu
        elif (( HASWELL == 1 ))
        then
                part=haswell
	else
                part=batch
	fi

        # Set the module type based on the input parameters
        if (( TWOTHOUSANDTWENTYONE == 1 ))
        then
                module=twothousand
        fi

        # Check if GROUPCUTOFF is enabled
        GROUPCUTOFF=$(head -5 "$1" | awk '{if($1=="#GROUPCUTOFF") print $2}')
        if (( GROUPCUTOFF == 1 ))
        then
                GROUP=1
        else
                GROUP=0
        fi

        # Create a new job submission file with the extracted parameters
        echo "#!/bin/bash -x" > SBATCHTMP
        echo "#SBATCH --nodes=1" >> SBATCHTMP
        echo "#SBATCH --ntasks=1" >> SBATCHTMP
        echo "#SBATCH --cpus-per-task=$NCPU" >> SBATCHTMP
        echo "#SBATCH --output=out.%j" >> SBATCHTMP
        echo "#SBATCH --error=err.%j" >> SBATCHTMP
        echo "#SBATCH --time="$TIME"" >> SBATCHTMP
        echo "#SBATCH --partition=$part" >> SBATCHTMP
  
        # Add GPU resources if needed
	if (( GPU == 1 ))
	then
		echo "#SBATCH --gres=gpu:1" >> SBATCHTMP
	fi	
        echo "source $tools/functions.sh" >> SBATCHTMP
	echo "module purge" >> SBATCHTMP
	if (( GPU == 1 ))
	then
		echo "module load gcc/7.3.1 cuda/9.2 fftw/gcc/single/sse/3.3.8 gromacs/nompi/gpu/gcc/single/2018.4" >> SBATCHTMP
	fi
	if (( HASWELL == 1 ))
	then
        	echo "module load intel/env/2018 intel/mpi/2018 fftw/intel/haswell/single/sse/3.3.8 gromacs/nompi/cpu/intel/haswell/single/2018.2" >> SBATCHTMP
	fi
	if (( TWOTHOUSANDTWENTYONE == 1 ))
	then
		if (( GPU == 1 ))
		then
			echo "module purge && module load gcc/7.3.1 fftw/gcc/single/sse/3.3.8 cuda/9.2 gromacs/nompi/gpu/gcc/single/2021.1" >> SBATCHTMP
		else
	        	echo "module purge && module load gcc/10.3.0 fftw/gcc/single/sse/3.3.9 gromacs/nompi/cpu/gcc/single/2021.1" >> SBATCHTMP
		fi
	else
        	echo "module load intel/env/2018 intel/mpi/2018 fftw/intel/single/sse/3.3.8 gromacs/nompi/cpu/intel/single/2018.4" >> SBATCHTMP
	fi
	if (( GROUP == 1 ))
	then
        	echo "MDRUN=\"gmx mdrun -ntmpi "$NCPU" -ntomp 1\"" >> SBATCHTMP
	else	
        	echo "MDRUN=\"gmx mdrun -ntomp "$NCPU" -ntmpi 1\"" >> SBATCHTMP
	fi
        echo " " >> SBATCHTMP
        echo "########################## " >> SBATCHTMP
        echo " " >> SBATCHTMP
        cat SBATCHTMP "$1" > SBTMP && mv SBTMP "$1" && rm SBATCHTMP
	chmod +x "$1"
}

# This function can be used for pull-off simulations in the z direction,
# as it makes the most efficient use of available cores by optimizing
# the domain decomposition and number of threads for the GROMACS mdrun command.
function update_mdrun_arg1_i_arg2_gro() {
    # Calculate rdd (real space cut-off) based on input argument 1 ($1)
    rdd=$(echo "2.5 + 0.1*$1" | bc -l)

    # Calculate rddscaled (scaled real space cut-off)
    rddscaled=$(echo "$rdd * 1.25" | bc -l)

    # Calculate xd (x-dimension decomposition) based on the last line of the input .gro file ($2)
    xd=$(echo "$(tail -1 $2 | awk -v "rdds=$rddscaled" '{print $1/rdds}') / 1" | bc)

    # Calculate yd (y-dimension decomposition) based on the last line of the input .gro file ($2)
    yd=$(echo "$(tail -1 $2 | awk -v "rdds=$rddscaled" '{print $2/rdds}') / 1" | bc)

    # Calculate the total number of threads
    threads=$(echo "$xd * $yd" | bc)

    # Calculate domain decomposition (dd) as a string containing xd, yd, and 1 (z-dimension)
    dd=$(echo "$xd $yd 1")

    # Check if the MDRUN variable starts with "gmx"
    bSNIFFA=$(echo $MDRUN | awk '{if($1=="gmx") print 1; else print 0}')

    # If MDRUN starts with "gmx", update the MDRUN variable with the calculated values
    if (( bSNIFFA ))
    then
        MDRUN=$(echo "gmx mdrun -ntmpi $threads -dd $dd" )
    # If MDRUN does not start with "gmx", update the MDRUN variable with the calculated values using gmx_mpi
    else
        MDRUN=$(echo "mpirun -np $threads gmx_mpi mdrun -dd $dd" )
    fi
}

