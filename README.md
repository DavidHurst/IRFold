I verify that I am the sole author of the programmes contained in this archive, except where explicitly stated to the contrary.
David Hurst, 12/08/2023

# IRFold

RNA Secondary Structure Prediction by optimising configurations of found inverted repeats via Constraint Programming

## Minimal Install
Assumes Conda installed 
1. Clone repository 
    ```
    $ git clone https://github.com/DavidHurst/IRFold 
    ```
2. Move to repository directory
    ```
    $ cd IRFold/
    ```
3. Create new environment from requirements file 
    ```
    $ conda env create -f requirements.yml
    ```
4. Activate environment 
    ```
    $ conda activate irfold
    ```
5. Run demo script 
    ```
    $ python demo.py
    ```

## Install for Experiments

1. Compile ViennaRNA
2. Compile GLPK
3. Compile IPknot  
5. Set location of thermodynamic parameters for RNAstructure:
    ```
    $ export DATAPATH=experiment_scripts/RNAstructure/data_tables/
    ```



