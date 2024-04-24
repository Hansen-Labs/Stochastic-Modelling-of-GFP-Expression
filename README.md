## GFP Expression Modeling
This repository contains the code used in the paper entitled "Decoupled translation and degradation enable noise-modulation by poly(A)-tail" by Grandi et al. (2024) to perform stochastic modeling of GFP expression from mRNAs with different poly(A)-tail lengths.

# Files
- model.py: This file contains the model and the rates needed for the simulation.
- simulate_independent.ipynb: Jupyter notebook used to execute the simulations. This notebook relies on model.py to perform the simulations.
  
# Output
The output of the simulations is stored in CSV files. These files contain the changes over time in the levels of RNA, premature GFP, and mature GFP for each simulated cell.
