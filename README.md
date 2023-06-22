# AIDPET

###  Automated Intrinsically Disordered Protein Equilibration Tool (for GROMACS) (AIDPET)

## Technology Highlights

This project effectively applies various domains, from molecular dynamics simulations to efficient automation pipelines. Some key highlights include:

 - Molecular Dynamics Simulations: We make proficient use of GROMACS and its included toolchain.
 - Intrinsically Disordered Proteins (IDPs): Statistical thermodynamic understanding of the unique properties of IDPs has allowed and informed the development of this specialized generator tailored for such proteins.
 - High-Performance Computing: Designed to run on HPC clusters using optimization techniques relevant to these computational environments.
 - Automation & Scripting: Leveraging Bash and Python scripting languages to create a robust and modular automation pipeline, streamlining the equilibration process and reducing manual intervention.
 - Strategic Approaches Employed: Effective implementation of error handling, amorphization techniques, and self-avoiding random walks in periodic boundary conditions, with attention to detail.

The codebase is an updated version of the code that has been used to publish ["How Does Gecko Keratin Stick to Hydrophilic and Hydrophobic Surfaces in the Presence and Absence of Water? An Atomistic Molecular Dynamics Investigation"](https://doi.org/10.1021/acsnano.2c08627). We explain the underlying concepts of molecular dynamics and computational biology in this readme but also as comments throughout the codebase. For more in-depth detail refer to the supporting information of the study mentioned above.

### If you use AIDPET in your work, please cite the following paper:

Materzok, T.; Canestraight, A.; Gorb, S.; Müller-Plathe, F. ["How Does Gecko Keratin Stick to Hydrophilic and Hydrophobic Surfaces in the Presence and Absence of Water? An Atomistic Molecular Dynamics Investigation"](https://doi.org/10.1021/acsnano.2c08627). ACS Nano 2022, 16 (11), 19261–19270.

## Overview

This repository provides AIDPET, a tool designed to automate a series of GROMACS molecular dynamics simulations for equilibrating complex protein systems and analyzing their mechanical properties under various conditions. The primary focus of this script is on gecko beta-associated-keratin proteins, aiming to automatically and robustly generate and analyze multiple gecko keratin protein sequences under different conditions.

## Key Features

AIDPET streamlines numerous essential steps in the process, such as:

 - Protein structure generation via multiple self-avoiding random walks in periodic boundary conditions
 - Energy minimization using soft-core potentials
 - High-temperature equilibration (amorphization) using soft-core potentials
 - High-temperature step-wise cooling
 - NPT relaxation
 - Disulfide bond formation
 - Water insertion
 - NPT equilibration with or without water
 - Stress-strain simulations considering different hydration and strain conditions

To ensure robustness and reliability, AIDPET is equipped with error-handling capabilities that promptly cancel all submitted jobs if any step encounters a failure. This comprehensive tool facilitates an automated pipeline of molecular dynamics simulations on an HPC Cluster and yields equilibrated IDP materials, including simulations of their stress-strain behavior. The equilibrated IDP material has a default cysteine cross-linking density of 33% and two systems will be created, with and without absorbed water.

## Intrinsically Disordered Proteins (IDPs)

AIDPET is particularly suited for studying intrinsically disordered proteins (IDPs), such as the core-box knockout sequences of the gecko beta keratin. Although proteins are generally not randomly ordered, IDPs exhibit a high degree of structural disorder and lack stable secondary and tertiary structures. This unique characteristic allows AIDPET to utilize amorphization techniques for generating IDP molecular dynamics models. Consequently, AIDPET may be applied to create and analyze other IDPs, extending its utility beyond gecko beta-associated keratin proteins.

## Dependencies

To use AIDPET, ensure the following dependencies are met:

- GROMACS 2018.4 and 2021.1 installed and available in the environment. 2021.1 is only used to speed-up some relaxation steps

- The python libraries are installed - Assuming the base is an Anaconda3 installation, this includes MDAnalysis and PeptideBuilder. 

- If 'slurm_present' in create_IDP_elastomer.sh is set to 'true', access to a Linux HPC cluster with SLURM is needed. If 'false', the script will directly run the simulation batch scripts with './job.sh'

- SAPGenPBC (https://github.com/TobiasMaterzok/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions) is present in 
  ~/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions
- AIDPET is in ~/AIDPET
- The scripts from AIDPET and SAPGenPBC link to the correct directories:
  - functions.sh
  - concatenate_pdbs.sh
  - create_amorph_ssbonds_lightH_modular.sh
  - create_gb36_disulfide_bondparameter.sh
  - functions.sh
  - get_frame_average_volume.sh
  - get_number_for_gmxenergy.sh
  - run_generators_and_concatenate_pdbs.sh
  - select_interactive_pdb2gmx_for_ssbonds_lightH.sh
  - select_interactive_pdb2gmx_lightH.sh
  - py_get_number_of_atomtype.py
  - py_numWater.py
  - SAPGenPBC.py
  
 - GROMACS .mdp files in the current working directory:
   - minim_first.mdp
   - minim.mdp
   - nvt_hot.mdp
   - mdp_nvt_cooldown.mdp
   - relax_angles_2.mdp
   - relax_angles_3_strong_w.mdp
   - minim_w_insert.mdp
   - enlarge_slow.mdp
   - enlarge.mdp

### Cloning the Repository and Setting Up

1. Clone the repository to your local machine, the scripts assumes this is done in the home directory:

```
cd ~
git clone https://github.com/TobiasMaterzok/AIDPET.git
cd AIDPET
```

2. Set the executable flag:

```
chmod +x *.sh
```

3. You need to also download the SAPGenPBC toolchain:

```
cd ~
git clone https://github.com/TobiasMaterzok/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions
cd PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions
chmod +x *.sh
```
again, assuming you are in the home directory.

4. Modify functions.sh and create_IDP_elastomer.sh to fit your HPC environment.

```
cd ~/AIDPET
```

5. Create your working directory:

```
mkdir -p /path/to/your/working/directory/
```

6. Copy the GROMACS .mdp files into your working directory:

```
cp gromacs_mdps/*.mdp /path/to/your/working/directory/
```

8. Make the scripts executable:

```
chmod +x ~/AIDPET/*.sh
chmod +x ~/PDB-Protein-Chain-Generator-in-Periodic-Boundary-Conditions/*.sh
```

Now you are all set to use AIDPET in your working directory.

## Usage

To run AIDPET, execute the following command in your working directory:

```
cd /path/to/your/working/directory/
~/AIDPET/create_IDP_elastomer.sh
```

If you run with slurm you can directly run AIDPET via:

```
cd /path/to/your/working/directory/
sbatch -J AIDPET ~/AIDPET/create_IDP_elastomer.sh
```

This will start the automation process for a series of GROMACS molecular dynamics simulations and yield a equilibrated protein systems at 1 bar and 300 K with 33% of the cysteines involved in cross-links in 3D periodic boundary conditions with 10 weight percent water and without water present.

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.    
