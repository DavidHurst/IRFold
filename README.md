I verify that I am the sole author of the programmes contained in this archive, except where explicitly stated to the contrary.
David Hurst, 12/08/2023

# Table of Contents

/data:
	Results of experiments and stored dataset (dataset omitted for sake of file size)
	
/experiment_scripts/experiment_1_*:
	Scripts which run and analyse experiment to show the additive free energy vs union free energy
/experiment_scripts/experiment_2_*:
	Scripts which run and analyse experiment to show that validating IR pairs also vaidates IR triplets and quadruplets
/experiment_scripts/experiment_3_*:
	Scripts which run and analyse experiment to show the error in the free energy additivity assumption
/experiment_scripts/experiment_4_*:
	Scripts which run and analyse experiment to show the performance of the solver for different IRFold variants on a test set of random sequences.
/experiment_scripts/experiment_5_*:
	Scripts which run and analyse experiment to benchmark the performance of IRFold against other models on the bpRNA-1m dataset
	
/irfold/irfold_base.py:
	Class defining the IRFoldBase model
/irfold/irfold_val1.py:
	Class defining the IRFoldVal1 model
/irfold/irfold_val2.py:
	Class defining the IRFoldVal2 model
/irfold/irfold_corx.py:
	Class defining the IRFoldCorX model
/irfold/irfold_cor2.py:
	Class used for development and testing of IRFoldCorX
/irfold/irfold_cor3.py:
	Class used for development and testing of IRFoldCorX
/irfold/IUPACpal:
	Binary exe of IUPACpal, used to find IRs
	
/test:
	Unit tests for package


# IRFold

RNA Secondary Structure Prediction by optimising configurations of found inverted repeats via Constraint Programming.
IRFold will only run on Linux due to RNAlib's Python package.

## Minimal Install
Assumes Conda installed 
1. Clone repository 
    ```
    git clone https://github.com/DavidHurst/IRFold 
    ```
2. Move to repository directory
    ```
    cd IRFold/
    ```
3. Create new environment from requirements file 
    ```
    conda env create -f requirements.yml
    ```
4. Activate environment 
    ```
    conda activate irfold
    ```
5. Run demo script 
    ```
    python demo.py
    ```

## Install for Experiments

1. Run minimal installation steps above
1. Compile and install ViennaRNA version 2.5.1 to /experiment_scripts/ViennaRNA-2.5.1 from https://www.tbi.univie.ac.at/RNA/
2. Compile and install GLPK version 5.0 to /experiment_scripts/ from https://www.gnu.org/software/glpk/
3. Compile and install IPknot to /experiment_scripts/ipknot from https://github.com/satoken/ipknot
5. Set location of thermodynamic parameters for RNAstructure:
    ```
    export DATAPATH=experiment_scripts/RNAstructure/data_tables/
    ```
6. Run experiment scripts of choice e.g.
    ```
    python experiment_scripts/experiment_1_investigate_ir_pair_fe_additivity.py
    ```

