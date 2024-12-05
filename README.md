PyCOREX – PyRosetta

These scripts are meant to be a tool and starting point for further exploration of the complex interface of protein structure prediction and design, through the lens
of physics, statistics, and thermodynamics.


COREX-Rosetta Structure Script Use

These notebooks require the use of two tools, the COREX algorithm and PyRosetta, a python accessible version of the Rosetta software suite. 

PyRosetta is developed out of the Gray lab at Johns Hopkins University and can be found with instructions for downloading to all operating systems here: https://www.pyrosetta.org

COREX is developed out of the Hilser lab at Johns Hopkins University and can be found with instructions for downloading to all operating systems here:


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



