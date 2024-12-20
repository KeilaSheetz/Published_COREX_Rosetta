PyCOREX – PyRosetta

These scripts are meant to be a tool and starting point for further exploration of the complex interface of protein structure prediction and design, through the lens
of physics, statistics, and thermodynamics.


COREX-Rosetta Structure Script Use

These notebooks require the use of two tools, the COREX algorithm and PyRosetta, a python accessible version of the Rosetta software suite. 

PyRosetta is developed out of the Gray lab at Johns Hopkins University and can be found with instructions for downloading to all operating systems here: https://www.pyrosetta.org

COREX is developed out of the Hilser lab at Johns Hopkins University.  COREX is written in C and must be compiled before use.  Although compiled code is available at this GitHub repository it is unlikely to work on your system (though feel free to try).  Stepwise intructions for compilation are the following, compiler gcc or equivalent is necessary.

1. Make a local subdirectory entited "COREX_Files" and download to there the contents of the GitHub COREX_Files subdirectory. 

2. Typing "make" within the shell will compile the Monte Carlo COREX ensemble generation program: EnsembleGeneratorMC.  (This program is used for medium size proteins, ~ 75 residues, or larger.)

3. Next, the Makefile must be edited by removing all instances of "MC" within the file.  Once this is done typing "make" again within the shell will compile the COREX full ensemble generation program used for small proteins: EnsembleGenerator.

4. Generated ensembles made from EnsembleGeneratorMC and EnsembleGenerator will be analyzed by the StructureFileReader programs, which must be compiled separately now in steps 4-7.  Type "gcc StructureFileReader_ThermoDescriptors_NativeMC.c -lm" within the shell.  The resulting "a.out" file should be renamed by typing "mv a.out a.out_nativeMC".

5. Similarly, type "gcc StructureFileReader_ThermoDescriptors_Native.c -lm" within the shell.  The resulting "a.out" file should be renamed by typing "mv a.out a.out_native".  

6. Similarly, type "gcc StructureFileReader_ThermoDescriptors_DenaturedMC.c -lm" within the shell.  The resulting "a.out" file should be renamed by typing "mv a.out a.out_denaturedMC".

7. Finally, type "gcc StructureFileReader_ThermoDescriptors_Denatured.c -lm" within the shell.  The resulting "a.out" file should be renamed by typing "mv a.out a.out_denatured".

8. Thus, the following compiled programs should exist in your COREX_Files subdirectory, as they will be called by the jupyter notebook: "EnsembleGeneratorMC", "EnsembleGenerator", "a.out_nativeMC", "a.out_native", "a.out_DenaturedMC", "a.out_Denatured".

9. These compiled programs may be used independently of the jupyter notebook, please see the README file for details.


SCORE NATIVE SEQUENCES ON IMPORTED STRUCTURES 

This script is for scoring experimentally solved structures with both COREX and Rosetta.

The only input for this script is a list of PDB codes. This import is made at the beginning of the script directly after the import cell - this is the only cell you are required to edit.

The output file will contain all the structures scored by Rosetta (in REU) 
and COREX (in Total LogOdds, also referred to as the Fold Recognition Score)


CREATE MUTANT SEQUENCES ON IMPORTED STRUCTURES 

This script is for threading any number of mutant sequences on to experimentally solved structures.

The inputs for this script are a list of PDB codes, and a list of names and mutant sequences to thread onto each of the structures. These imports are made at the beginning of the script directly after the import cell. The only cell you need to edit is this cell.

IMPORTANT: All PDBs and sequences must be exactly the same length for this version of the script to work. 

The output file will contain all the mutant sequences scored on all desired folds from Rosetta (in REU) and COREX (in Total LogOdds, also referred to as the Fold Recognition Score)



