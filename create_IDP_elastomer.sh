#!/bin/sh -x
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH --time=4-00:00:00
#SBATCH --partition=batch # Please modify according to your partitions

#
# Equilibrated IDP (Protein) System Generator
# 
# Author: Tobias Materzok https://github.com/TobiasMaterzok/
#
# This script automates a series of GROMACS molecular dynamics simulations for a the equilibration
# of a complex protein system and analysis of their mechanical properties under various conditions.
# The system under study is a gecko beta-associated-keratin protein, and 
# the simulations are aimed at generating and analyzing various conformations of the protein under 
# different conditions. The script performs several steps, including protein structure generation, 
# energy minimization, temperature equilibration, cooling from high temperatures, NPT relaxation,
# disulfide bond formation, water insertion, NPT equilibration with and without water inserted, 
# and stress-strain simulations under varying hydration and strain conditions.
# 
# The script also includes error handling code to cancel all submitted jobs if any step fails.
#
# Usage: ./create_IDP_elastomer.sh
#
# Dependencies:
# - GROMACS 2018.4 and 2021.1 installed and available in the environment
# - SAPGenPBC (https://github.com/TobiasMaterzok/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions) is present in 
#   ~/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions
# - AIDPET is in ~/AIDPET
# - The scripts from AIDPET and SAPGenPBC link to the correct directories:
# - - concatenate_pdbs.sh
# - - create_amorph_ssbonds_lightH_modular.sh
# - - create_gb36_disulfide_bondparameter.sh
# - - functions.sh
# - - get_frame_average_volume.sh
# - - get_number_for_gmxenergy.sh
# - - run_generators_and_concatenate_pdbs.sh
# - - select_interactive_pdb2gmx_for_ssbonds_lightH.sh
# - - select_interactive_pdb2gmx_lightH.sh
# - - py_get_number_of_atomtype.py
# - - py_numWater.py
# - - SAPGenPBC.py
# - GROMACS .mdp files in the current working directory:
# - - minim_first.mdp
# - - minim.mdp
# - - nvt_hot.mdp
# - - mdp_nvt_cooldown.mdp
# - - relax_angles_2.mdp
# - - relax_angles_3_strong_w.mdp
# - - minim_w_insert.mdp
# - - enlarge_slow.mdp
# - - enlarge.mdp
# - If 'slurm_present' is set to 'true', access to a Linux HPC cluster with SLURM installed and 
# configured is needed. If 'false', the script will directly run the simulation batch scripts with './job.sh'
# Input: PDB structure of the protein, box size and density specified in the script; 
#

# Adjust the 'corenum' variable and the '--cpus-per-task' parameter according to the resources available on the HPC cluster
corenum=$SLURM_CPUS_PER_TASK
MDRUN="gmx mdrun -ntomp $corenum -ntmpi 1"

source ~/anaconda3/etc/profile.d/conda.sh
# conda create --name role_of_seq
# conda install -c conda-forge mdanalysis
# pip install biopython
# pip install PeptideBuilder
conda activate role_of_seq

# Specify path to AIDPET
AIDPET=~/AIDPET
# Specify path to SAPGenPBC 
SAPGenPBC=~/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions
source $AIDPET/functions.sh

load_gromacs_2018_4

use_run_generators=true
# slurm_present: Set this variable to 'true' if on an HPC cluster with SLURM. 
# If 'false', the script will directly run the simulation batch scripts with './job.sh' instead of using 'sbatch job.sh'.
slurm_present=true

parentdir=$(pwd)
for L in lseq_test_70_70_70_dens12_sample_4
do
	jobname="$L"
	mkdir "$parentdir"/"$L"
	cd "$parentdir"/"$L"
	simpath=$(pwd)

	# The if-else statement in this automation script checks for the success of the previous step
	# by verifying the existence of a specific GROMACS output file (gro). If the file exists, it 
	# assumes success and proceeds with the next simulation. If not, it cancels all submitted jobs 
	# in the current $simpath directory, using job IDs stored in the 'X' variable, indicating the 
	# previous step's failure.
	if [[ "$use_run_generators" = "true" ]]
	then
		$SAPGenPBC/run_generators_and_concatenate_pdbs.sh $corenum MSCCPPSCATPSCPKPCCSPCCSPCGYPTGGLGSLGCCPCPCGPSSCCGSSTSARCLGITSGASVSCINQIPASCEPLRVGGYTACGGCPPCGRIC 70 70 70 1.2
	else
		cp ../*.pdb .
	fi

	numchains=$(ls -l Chain*.pdb | wc -l)
	cp ../*mdp .

	# Automatically generate GROMACS .gro and .top file for the current protein system using the GROMOS 54A7 force field
	$AIDPET/select_interactive_pdb2gmx_lightH.sh "GA9rnd_gromos.conf" "GA9rnd_top_gromos.top"
	cp GA9rnd_gromos.conf GA9rnd_gromos.conf_backup
	cp GA9rnd_top_gromos.top GA9rnd_top_gromos.top_backup

	cp Chain1.pdb Chainreference.pdb
	mkdir Chains
	mv Chain?.pdb Chains/
	mv Chain??.pdb Chains/
	mv Chain???.pdb Chains/


	# Step 1: Energy Minimization
	# We create and configure the job script 'job_firstpart_minim.sh' for the minimization process.
	# The purpose of energy minimization is to remove any steric clashes or inappropriate geometry in the starting
	# structure. This is done here, using the Steepest Descent algorithm to minimize the potential energy of the system.
	# We perform energy minimization with both H-bonds and all-bonds constraints. The number of steps for the first
	# minimization is set to 300. We then loop over lambda values from 5 to 20 with a step of 1, updating the
	# 'init_lambda_state' parameter in the 'minim_first.mdp' file for each iteration. This ensures a gradual decrease
	# in potential energy due to the soft-core potentials, reducing the chance of local minima trapping.
	# We submit the job script to the scheduler using 'sbatch'
	# and start creating a list of dependencies using the naming scheme "$jobname"_"$runnum"
	runnum=1
	currname=job_firstpart_minim.sh
	rm $currname
	batch_template $currname $corenum "GROUPCUTOFF 0" "2021 1"
	change_sbatch $currname
	change_mdp constraints h-bonds minim_first.mdp #change_mdp constraints none minim_first.mdp  "NONE for samples 1-3" now h-bond for sample 4
	change_mdp constraints all-bonds minim.mdp
	change_mdp nsteps 500 minim_first.mdp
	echo "cp GA9rnd_gromos.conf.gro GA9rnd_gromos_em_first.gro" >> $currname
	echo "for i in \$(seq 5 1 20)" >> $currname
	echo "do" >> $currname
	echo "change_mdp init_lambda_state \$i minim_first.mdp"  >> $currname
	echo "gmx grompp -f minim_first.mdp -c GA9rnd_gromos_em_first.gro -p GA9rnd_top_gromos.top -o GA9rnd_gromos_em_first.tpr -maxwarn 1" >> $currname
	echo "\$MDRUN -v -deffnm GA9rnd_gromos_em_first" >> $currname
	echo "rm \#*" >> $currname
	echo "done" >> $currname
	echo "gmx grompp -f minim.mdp -c GA9rnd_gromos_em_first -p GA9rnd_top_gromos.top -o GA9rnd_gromos_em.tpr -maxwarn 1" >> $currname
	echo "\$MDRUN -v -deffnm GA9rnd_gromos_em" >> $currname
	
	if [[ "$slurm_present" == "true" ]]; then 
		sbatch -J "$jobname"_"$runnum" "$currname"
	else
		chmod +x $currname
		./"$currname"
	fi

	# Step 2: Temperature Equilibration (Hot)
	# The hot temperature equilibration is performed to quickly relax any unrealistic frustrations on the picosecond timescale,
	# allowing the chains to lose any correlations to their initial conformations. This is due to the soft-core potentials and
	# the high annealing temperature of 1300 K (RT = 11 kJ/mol). The decay of the autocorrelation function of the end-to-end vector
	# and the radius of gyration is reduced to approximately 25 ps. We use a range of lambda values from 0.7 to 1.0 in steps of 0.1,
	# performing four runs of 100 ps each with a timestep of 0.4 fs. This "slow push-off" procedure is known to reduce
	# perturbations in local chain conformations. We use a relative dielectric constant of 80 to simulate the electrostatic
	# screening effect of water, and a Berendsen thermostat with a coupling time constant of 0.1 ps in the NVT ensemble.
	mkdir nvt_hot
	cd nvt_hot
		runprev=$runnum
		((runnum++))
		currname=job_secondpart_hot.sh
		rm $currname
		batch_template "$currname" $corenum "GROUPCUTOFF 0" "2021 1"
		change_sbatch "$currname"
		cp ../nvt_hot.mdp .
		echo "cp ../GA9rnd_gromos_em.gro GA9rnd_gromos_hot.gro" >> $currname
		echo "for i in \$(seq 14 2 20)" >> $currname
		echo "do" >> $currname
		echo "change_mdp init_lambda_state \$i nvt_hot.mdp"  >> $currname
		echo "time=0.1" >> $currname
		echo "dt=0.0004" >> $currname
		echo "steps=\$(ns_to_steps \$time \$dt)" >> $currname
		echo "change_mdp nsteps \$steps nvt_hot.mdp" >> $currname
		echo "change_mdp dt \$dt nvt_hot.mdp" >> $currname
		echo "gmx grompp -f nvt_hot.mdp -c GA9rnd_gromos_hot.gro -p ../GA9rnd_top_gromos.top -o GA9rnd_gromos_hot.tpr -maxwarn 1" >> $currname
		echo "\$MDRUN -v -deffnm GA9rnd_gromos_hot" >> $currname
		echo "cp GA9rnd_gromos_hot.gro ../run_nvt_cooldown_confout_NPT.gro" >> $currname
		echo "done" >> $currname
		echo "X=\$(cat ../IDs_FOR_CANCELING_IF_ERROR)" >> $currname
		echo "if [ -f GA9rnd_gromos_hot.gro ]; then echo YO; else scancel \"\$X\"; echo \"FAILED at Energy minimization\" >> ../ERROR; fi" >> $currname
		if [[ "$slurm_present" == "true" ]]; then 
			sbatch -J "$jobname"_"$runnum" --dependency=$(squeue --noheader --format %i --name "$jobname"_"$runprev") "$currname"
		else
			chmod +x $currname
			./"$currname"
		fi
	cd $simpath

	# Step 3: Cooldown from 1300K to 300K
	# In this step, we gradually decrease the temperature from 1300K to 300K in increments of -50K. This cooldown process
	# allows the protein chain to explore more energetically favorable conformations and helps in achieving better
	# sampling of the conformational space. At each temperature level, we perform NVT simulations with different time
	# steps and number of steps depending on the temperature, using Berendsen thermostat and, when necessary,
	# Berendsen pressure coupling. This helps in maintaining stable temperature and pressure while reducing the
	# simulation time, resulting in a more efficient exploration of the conformational space.
	mkdir nvt_cooldown_NPT
	cd nvt_cooldown_NPT
		for i in $(seq 1300 -50 300)
		do
		mkdir nvt_cooldown_"$i"
		cd nvt_cooldown_"$i"
			steps=$(ns_to_steps 0.25 0.0005)
			cp ../../mdp_nvt_cooldown.mdp mdp_nvt_cooldown_"$i".mdp
			change_mdp constraints h-bonds mdp_nvt_cooldown_"$i".mdp
			change_mdp nsteps $steps mdp_nvt_cooldown_"$i".mdp
			change_mdp dt 0.0005 mdp_nvt_cooldown_"$i".mdp
			change_mdp ref_t $i mdp_nvt_cooldown_"$i".mdp

			change_mdp epsilon-r $(echo "$i/100*8-23" | bc) mdp_nvt_cooldown_"$i".mdp
			change_mdp gen-vel yes mdp_nvt_cooldown_"$i".mdp
			change_mdp gen-temp $i mdp_nvt_cooldown_"$i".mdp

			# Here we adjust the GROMACS simulation parameters based on the current iteration number, 
			# which corresponds to the current temperature in the cooling protocol. We increase the 
			# timestep and, thus, can reduce the number of simulation steps as the system cools down 
			# and the time evolution slows down, optimizing simulation performance while maintaining accuracy. 
			# We also add a barostat at lower temperatures. Overall, these adjustments to the simulation parameters 
			# are to speedup each stage of the cooling protocol, resulting in faster and robust simulations
			if (( $i < 800 )); then
				steps=$(ns_to_steps 0.25 0.0006)
				change_mdp nsteps $steps mdp_nvt_cooldown_"$i".mdp
				change_mdp dt 0.0006 mdp_nvt_cooldown_"$i".mdp
			fi
			if (( $i < 700 )); then
				steps=$(ns_to_steps 0.25 0.001)
				change_mdp nsteps $steps mdp_nvt_cooldown_"$i".mdp
				change_mdp dt 0.001 mdp_nvt_cooldown_"$i".mdp
			fi
			if (( $i < 400 )); then
				steps=$(ns_to_steps 2.5 0.002)
				change_mdp nsteps $steps mdp_nvt_cooldown_"$i".mdp
				change_mdp dt 0.002 mdp_nvt_cooldown_"$i".mdp
				change_mdp pcoupl berendsen mdp_nvt_cooldown_"$i".mdp
				change_mdp tau_p 0.5 mdp_nvt_cooldown_"$i".mdp
			fi
			if [ $i == 300 ]; then 
				steps=$(ns_to_steps 5 0.002)
				change_mdp nsteps $steps mdp_nvt_cooldown_"$i".mdp
				change_mdp dt 0.002 mdp_nvt_cooldown_"$i".mdp
			fi

			runprev=$runnum
			((runnum++))
			currname=job_nvt_cooldown_"$i".sh
			rm $currname
			batch_template $currname $corenum "GROUPCUTOFF 0"
			change_sbatch $currname
			echo "gmx grompp -f mdp_nvt_cooldown_"$i".mdp -c ../../run_nvt_cooldown_confout_NPT.gro -p ../../GA9rnd_top_gromos.top -o run_nvt_cooldown_"$i"_NPT.tpr -maxwarn 1" >> $currname
			echo "\$MDRUN -deffnm run_nvt_cooldown_"$i"_NPT -v" >> $currname
			echo "cp run_nvt_cooldown_"$i"_NPT.gro ../../run_nvt_cooldown_confout_NPT.gro" >> $currname
			
			if [[ "$slurm_present" == "true" ]]; then 
				sbatch -J "$jobname"_"$runnum" --dependency=$(squeue --noheader --format %i --name "$jobname"_"$runprev") "$currname"
			else
				chmod +x $currname
				./"$currname"
			fi
		cd ..
		done
	cd $simpath

	# Step 4: Relaxation of Backbone Angles under NPT Conditions
	# The purpose of this step is to mainly relax the long-range backbone angles while converging to the 
	# equilibrium inter protein-protein distances. This is done by running a 10 ns NPT simulation with a 
	# 2 fs time step using all-bonds constraints, a pressure coupling constant (tau_p) of 0.5 ps, and a 
	# temperature coupling constant (tau_t) of 0.1 ps. The simulation is performed in the 'relax_angles_NPT' 
	# directory. This step allows for further exploration of the conformational space while preparing the 
	# system for disulfide bond formation in the next step.
	mkdir relax_angles_NPT
	cd relax_angles_NPT
		runprev=$runnum
		((runnum++))
		currname=job_relax_angles_tm.sh
		rm $currname
		batch_template $currname $corenum "GROUPCUTOFF 0"
		change_sbatch $currname
		cp ../relax_angles_2.mdp relax_angles.mdp
		steps=$(ns_to_steps 10 0.002)
		change_mdp nsteps $steps relax_angles.mdp
		change_mdp dt 0.002 relax_angles.mdp
		change_mdp constraints all-bonds relax_angles.mdp
		change_mdp tau_p 0.5 relax_angles.mdp
		change_mdp tau_t 0.1 relax_angles.mdp
		echo "gmx grompp -f relax_angles.mdp -c ../run_nvt_cooldown_confout_NPT.gro -p ../GA9rnd_top_gromos.top -o GA9rnd_gromos_relaxangles_2_NPT.tpr -maxwarn 1" >> $currname
		echo "\$MDRUN -v -deffnm GA9rnd_gromos_relaxangles_2_NPT" >> $currname
		echo "cp GA9rnd_gromos_relaxangles_2_NPT.gro ../" >> $currname
		
		if [[ "$slurm_present" == "true" ]]; then 
			sbatch -J "$jobname"_"$runnum" --dependency=$(squeue --noheader --format %i --name "$jobname"_"$runprev") "$currname"
		else
			chmod +x $currname
			./"$currname"
		fi
		
	cd $simpath
	atomid=$(tail -3 Chainreference.pdb | head -1 | awk '{print $2}')
	atomtype=$(tail -3 Chainreference.pdb | head -1 | awk '{print $3}')
	resnumber=$(tail -3 Chainreference.pdb | head -1 | awk '{print $6}')

	runprev=$runnum
	runbeforess=$runprev

	# Step 5: Disulfide Bond Formation
	# In this step, we aim to create disulfide bonds between cysteine residues. This is important for achieving the target
	# disulfide bond density in the protein, which in the case of gecko beta keratin is approximately 7.5% (33% of all cysteines
	# being involved in a disulfide bridge). The "create_amorph_ssbonds_lightH_modular.sh" script is utilized to
	# generate the number of disulfide bonds required by iteratively increasing a cutoff distance for disulfide bond
	# formation. The process continues until the target number of disulfide bonds is reached or the maximum number of
	# iterations (60) is exceeded. The NPT simulation system prepared in the previous step is used as input for this
	# disulfide bond formation process, ensuring that the protein conformations are already closer to their equilibrium
	# protein-protein distances.
	for sdens in 0.33
	do
		((runnum++))
		rm -r $simpath/ss_"$sdens"
		mkdir $simpath/ss_"$sdens"
		cd $simpath/ss_"$sdens"
		cp ../*mdp .

		currname=job_sshbonds_NPT.sh
		batch_template $currname $corenum "GROUPCUTOFF 0"
		change_sbatch $currname
		echo "cp ../relax_angles_NPT/GA9rnd_gromos_relaxangles_2_NPT.* ." >> $currname
		echo "cp ../*top ." >> $currname
		echo "cp ../*ndx ." >> $currname
		echo "X=\$(cat IDs_FOR_CANCELING_IF_ERROR)" >> $currname
		echo "if [ -f GA9rnd_gromos_relaxangles_2_NPT.gro ]; then" >> $currname
		echo "cp relax_angles_NPT/GA9rnd_gromos_relaxangles_2_NPT.* ." >> $currname
		echo "cp GA9rnd_gromos_relaxangles_2_NPT.gro GA9rnd_gromos_relaxed_NPT.gro" >> $currname
		
		# Copy and modify the "create_amorph_ssbonds_lightH_modular.sh" script to create disulfide bonds in the NPT simulation
		cp $AIDPET/create_amorph_ssbonds_lightH_modular.sh create_amorph_ssbonds_NPT.sh
		chmod +x create_amorph_ssbonds_NPT.sh
		change_sbatch create_amorph_ssbonds_NPT.sh
		
		# Execute the disulfide bonds creation script with required arguments
		echo "./create_amorph_ssbonds_NPT.sh GA9rnd_gromos_relaxangles_2_NPT $resnumber $sdens" >> $currname
		
		# Create and populate "index.ndx" file with disulfide bond information in the "ssbond" directory
		echo "echo \"[ Link ]\" > ssbond/index.ndx" >> $currname
		echo "cat ssbond/sel_SS_tmp.dat | awk '{printf(\"%s \",substr(\$2,5,6))}' >> ssbond/index.ndx" >> $currname
		echo "cat ssbond/sel_SS_tmp.dat | awk '{printf(\"%s \",substr(\$5,5,6))}' >> ssbond/index.ndx" >> $currname
		echo "echo \" \" >> ssbond/index.ndx" >> $currname
		
		echo "else scancel \"\$X\"; echo \"FAILED before ssbond creation\" >> ERROR; fi" >> $currname

		if [[ "$slurm_present" == "true" ]]; then 
			sbatch -J "$jobname"_"$runnum"_"$sdens" --dependency=$(squeue --noheader --format %i --name "$jobname"_"$runbeforess") "$currname" 
		else
			chmod +x $currname
			./"$currname"
		fi
		
		# Step 6: NPT Relaxation of Crosslinked Amorphous Protein Material
		# The purpose of this step is to quickly relax the crosslinked amorphous protein material under NPT conditions, allowing for
		# the protein structure. This step is performed after disulfide bond formation to ensure that the newly created bonds are 
		# properly incorporated into the system. The simulation is carried out for 10 ns with a time step of 4 fs using all-bonds
		# constraints, a pressure coupling constant (tau_p) of 0.5 ps, and a temperature coupling constant (tau_t) of 0.1 ps. 
		# Prior to the NPT relaxation, the system undergoes energy minimization with both H-bonds and all-bonds constraints, using a
		# loop over (soft-core) lambda values from 5 to 20. This, again, ensures a gradual decrease in potential energy and reduces 
		# the chance of local minima trapping. The NPT relaxation is performed in the 'npt_relax_NPT' directory.
		mkdir npt_relax_NPT
		cd npt_relax_NPT
			runprev=$runnum
			((runnum++))
			currname=job_npt_relax.sh
			batch_template $currname $corenum "GROUPCUTOFF 0"
			change_sbatch $currname
			cp ../relax_angles_2.mdp npt_relax_after_ss.mdp
			steps=$(ns_to_steps 10 0.002)
			change_mdp nsteps $steps npt_relax_after_ss.mdp
			change_mdp dt 0.002 npt_relax_after_ss.mdp
			change_mdp constraints all-bonds npt_relax_after_ss.mdp
			change_mdp tau_p 0.5 npt_relax_after_ss.mdp
			change_mdp tau_t 0.1 npt_relax_after_ss.mdp
			cp ../minim.mdp .
			echo "periodic-molecules       = yes" >> minim.mdp
			change_mdp constraints all-bonds minim.mdp
			cp ../minim_first.mdp .
			change_mdp constraints h-bonds minim_first.mdp
			echo "periodic-molecules       = yes" >> minim_first.mdp
			change_mdp coulombtype PME minim_first.mdp
			change_mdp "epsilon-r" 1 minim_first.mdp
			change_mdp nsteps 3000 minim_first.mdp
			change_mdp bonded_lambdas "0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00" minim_first.mdp
			echo "cp ../ssbond/GA9rnd_gromos_ss.top ../GA9rnd_gromos_ss_NPT.top" >> $currname
			echo "for i in \$(seq 5 1 20)" >> $currname
			echo "do" >> $currname
			echo "change_mdp init_lambda_state \$i minim_first.mdp"  >> $currname
			echo "gmx grompp -f minim_first.mdp -c ../ssbond/GA9rnd_gromos_ss.gro -p ../GA9rnd_gromos_ss_NPT.top -o GA9rnd_gromos_ss_em_first.tpr -maxwarn 1" >> $currname
			echo "\$MDRUN -v -deffnm GA9rnd_gromos_ss_em_first" >> $currname
			echo "done" >> $currname
			echo "rm \#*" >> $currname
			echo "gmx grompp -f minim.mdp -c GA9rnd_gromos_ss_em_first.gro -p ../GA9rnd_gromos_ss_NPT.top -o GA9rnd_gromos_ss_em.tpr -maxwarn 1" >> $currname
			echo "\$MDRUN -v -deffnm GA9rnd_gromos_ss_em" >> $currname
			echo "gmx grompp -f npt_relax_after_ss.mdp -c GA9rnd_gromos_ss_em.gro -p ../GA9rnd_gromos_ss_NPT.top -o GA9rnd_npt_relax.tpr -maxwarn 1" >> $currname
			echo "\$MDRUN -v -deffnm GA9rnd_npt_relax" >> $currname
			
			if [[ "$slurm_present" == "true" ]]; then 
				sbatch -J "$jobname"_"$runnum"_"$sdens" --dependency=$(squeue --noheader --format %i --name "$jobname"_"$runprev"_"$sdens") "$currname"
			else
				chmod +x $currname
				./"$currname"
			fi
			
		cd $simpath/ss_"$sdens"
	
	
		runprev=$runnum
		runprev_waterrelax=$runprev
	
		# Step 7: NPT relaxation with water insertion
		# 1. Set the number of steps and simulation parameters for the next relaxation phase to 75 ns.
		# 2. Add periodic-molecules constraint to the minim_w_insert.mdp file.
		# 3. Create a directory for NPT relaxation with water inserted.
		# 4. Loop over water content (wt) values of 0.0 and 0.1 (without and with water).
		# 5. Prepare the input files and perform energy minimization with water inserted.
		# 6. Run NPT relaxation simulations.
		# 7. Calculate the average volume of the simulation box over the last 17% of the trajectors using get_frame_average_volume.sh. 
		#    This script selects a representative configuration with a volume closest to the average volume, ensuring that the 
		#    system has experienced a sufficiently long equilibration time and is representative of the equilibrium state.
		# 8. Submit the job with a dependency on the previous water relaxation step.
		steps=$(ns_to_steps 75 0.002)
		change_mdp nsteps $steps relax_angles_3_strong_w.mdp
		change_mdp dt 0.002 relax_angles_3_strong_w.mdp
		change_mdp tau_p 0.5 relax_angles_3_strong_w.mdp
		change_mdp tau_t 0.1 relax_angles_3_strong_w.mdp
		echo "periodic-molecules       = yes" >> minim_w_insert.mdp
		change_mdp constraints all-bonds minim_w_insert.mdp
		change_mdp constraints all-bonds relax_angles_3_strong_w.mdp
		mkdir npt_relax_water_NPT
		cd npt_relax_water_NPT
			for wt in 0.0 0.1
			do
			((runnum++))
			runprev_waterrun=$runnum
			mkdir $simpath/ss_"$sdens"/npt_relax_water_NPT/npt_relax_w_"$wt"
			cd $simpath/ss_"$sdens"/npt_relax_water_NPT/npt_relax_w_"$wt"
	
			currname=job_npt_relax_water.sh
			batch_template $currname $corenum "GROUPCUTOFF 0"
			change_sbatch $currname
	
			echo "cp ../../GA9rnd_gromos_ss_NPT.top GA9rnd_gromos_ss_w_"$wt"_NPT.top" >> $currname
			echo "cp ../../npt_relax_NPT/GA9rnd_npt_relax.gro GA9rnd_gromos_ss_w.gro" >> $currname
			echo "numW=\$(python $AIDPET/py_numWater.py ../../npt_relax_NPT/GA9rnd_npt_relax.gro $wt)" >> $currname
	
			echo "if (( \$numW > 0 )); then" >> $currname
			echo "gmx insert-molecules -f ../../npt_relax_NPT/GA9rnd_npt_relax.gro -o GA9rnd_gromos_ss_w.gro -ci ~/spc.gro -nmol "\$numW" -try 30000 -scale 0.70" >> $currname
			echo "echo "SOL \$numW" >> GA9rnd_gromos_ss_w_"$wt"_NPT.top" >> $currname
			echo "fi" >> $currname
	
			echo "gmx grompp -f ../../minim_w_insert.mdp -c GA9rnd_gromos_ss_w.gro -p GA9rnd_gromos_ss_w_"$wt"_NPT.top -o GA9rnd_gromos_ss_em_w.tpr" >> $currname
			echo "\$MDRUN -v -deffnm GA9rnd_gromos_ss_em_w" >> $currname
			echo "gmx grompp -f ../../relax_angles_3_strong_w.mdp -c GA9rnd_gromos_ss_em_w.gro -p GA9rnd_gromos_ss_w_"$wt"_NPT.top -o GA9rnd_npt_relax_strong_w_"$wt".tpr" >> $currname
			echo "\$MDRUN -v -deffnm GA9rnd_npt_relax_strong_w_"$wt"" >> $currname
			
			nstcalcenergy=$(awk '{if($1=="nstcalcenergy"){print $3}}' ../../relax_angles_3_strong_w.mdp)
			nstxoutcompressed=$(awk '{if($1=="nstxout-compressed"){print $3}}' ../../relax_angles_3_strong_w.mdp)
			arg2=$(python -c "print($nstxoutcompressed//$nstcalcenergy)")
			echo "$AIDPET/get_frame_average_volume.sh GA9rnd_npt_relax_strong_w_"$wt" $arg2" >> $currname
			
			if [[ "$slurm_present" == "true" ]]; then 
				sbatch -J "$jobname"_"$runnum"_"$sdens" --dependency=$(squeue --noheader --format %i --name "$jobname"_"$runprev_waterrelax"_"$sdens") "$currname"
			else
				chmod +x $currname
				./"$currname"
			fi
			
			cd $simpath/ss_"$sdens"/npt_relax_water_NPT
			
			runprev=$runprev_waterrun
			((runnum++))


			# Step 8: Stress-strain simulations
			# In this step, we perform stress-strain simulations for different water contents (wt) and strain values. 
			# The purpose of these simulations is to analyze the mechanical properties of the system under varying conditions 
			# of hydration and applied strain.
			# 1. For each water content (wt), a directory named "stressstrain_$wt" is created to store the respective simulation files.
			# 2. The MDP files "enlarge_slow.mdp" and "enlarge.mdp" are copied from the parent directory and modified 
			#    by setting constraints, tau_p, and tau_t parameters.
			# 3. The code iterates over strain values (0 to 2) with increments of 1, and for each strain value, 
			#    a subdirectory named "strain_simulation_$i" is created.
			# 4. If the strain value is greater than 0, a non-equilibrium molecular dynamics (NEMD) simulation using 
			#    the slow straining protocol ("enlarge_slow.mdp") is executed.
			# 5. An equilibrium molecular dynamics (EMD) simulation using the "enlarge.mdp" file is executed for all strain values.
			# 6. The EMD simulation takes the output configuration and topology files from the previous NPT relaxation step 
			#    with water inserted, and scales the system accordingly based on the calculated strain value.
			# 7. Both the NEMD and EMD simulation jobs are submitted with a dependency on the previous simulation, 
			#    ensuring the correct execution order.
			#
			# In summary, this step performs stress-strain simulations for different hydration and strain conditions, 
			# allowing for the analysis of the mechanical properties of the system under various conditions.
			mkdir stressstrain_"$wt"
			cd stressstrain_"$wt"
				cp $simpath/../enlarge_slow.mdp .
				cp $simpath/../enlarge.mdp .
				change_mdp constraints all-bonds enlarge.mdp
				change_mdp constraints all-bonds enlarge_slow.mdp
				change_mdp tau_p 0.5 enlarge.mdp
				change_mdp tau_t 0.1 enlarge.mdp
				change_mdp tau_p 0.5 enlarge_slow.mdp
				change_mdp tau_t 0.1 enlarge_slow.mdp

				for i in $(seq 0 1 2)
				do
					((runnum++))
					mkdir strain_simulation_"$i"
					cd strain_simulation_"$i"
	
					if [[ $i -gt 0 ]];
					then
						currname=job_stressstrain_nemd.sh
						batch_template $currname $corenum "GROUPCUTOFF 0"
						change_sbatch $currname
						echo "cp ../../npt_relax_w_"$wt"/confout_avgvol.gro ." >> $currname
						echo "cp ../../npt_relax_w_"$wt"/GA9rnd_gromos_ss_w_"$wt"_NPT.top ." >> $currname
						echo "gmx grompp -f ../enlarge_slow.mdp -c confout_avgvol.gro -p GA9rnd_gromos_ss_w_"$wt"_NPT.top -o run_strain.tpr -maxwarn 1" >> $currname
						echo "\$MDRUN -v -deffnm run_strain" >> $currname
						
						if [[ "$slurm_present" == "true" ]]; then 
							sbatch -J "$jobname"_"$runnum"_"$sdens" --dependency=$(squeue --noheader --format %i --name "$jobname"_"$runprev"_"$sdens") "$currname"
						else
							chmod +x $currname
							./"$currname"
						fi
						
						((runnum++))
					fi
					currname=job_stressstrain_eq.sh
					batch_template $currname $corenum "GROUPCUTOFF 0"
					change_sbatch $currname
					echo "cp ../../npt_relax_w_"$wt"/confout_avgvol.gro ." >> $currname
					echo "cp ../../npt_relax_w_"$wt"/GA9rnd_gromos_ss_w_"$wt"_NPT.top ." >> $currname
					echo "strainnum=$i" >> $currname
					echo "strain=\$(echo \"1+\$strainnum/100.0\" | bc -l)" >> $currname
					echo "gmx editconf -f confout_avgvol.gro -scale 1.0 1.0 \"\$strain\" -o run_strain_\"\$strainnum\".gro" >> $currname
					echo "gmx grompp -f ../enlarge.mdp -c run_strain_\"\$strainnum\".gro -p GA9rnd_gromos_ss_w_"$wt"_NPT.top -o run_strain_"$i".tpr -maxwarn 1" >> $currname
					echo "\$MDRUN -v -deffnm run_strain_"$i"" >> $currname
					
					if [[ "$slurm_present" == "true" ]]; then 
						sbatch -J "$jobname"_"$runnum"_"$sdens" --dependency=$(squeue --noheader --format %i --name "$jobname"_"$runprev"_"$sdens") "$currname"
					else
						chmod +x $currname
						./"$currname"
					fi
										
					cd ..
				done
			done
			cd $simpath/ss_"$sdens"/npt_relax_water_NPT
			runprev=$runnum
			((runnum++))
		cd $simpath/ss_"$sdens"
	done	

	# The loop iterates over the range of submitted job numbers (from 1 to 'runnum'), 
	# extracting the respective job IDs using the 'squeue' command with the '--noheader' and '--format' options.
	# The extracted job IDs are piped to AWK to format them as a space-separated list of IDs. 
	# Finally, the list of job IDs is written to a file named "IDs_FOR_CANCELING_IF_ERROR", 
	# which can be used to quickly cancel all submitted jobs in case of an error or if the user 
	# needs to stop the simulations.
	runprev=$runnum
	((runnum++)) 
	if [[ "$slurm_present" == "true" ]]; then
		for K in $(seq 1 1 $runnum)
		do
			ID=$(squeue --noheader --format %i --name "$jobname"_"$K")
			printf "$ID \n"
		done | awk '{printf("%s ",$1)}' > IDs_FOR_CANCELING_IF_ERROR
	fi
cd $parentdir
done
