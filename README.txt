This repository includes the distribution data and MATLAB code used to analyze this data, for Norris, et al. "Bacterial chemotaxis to saccharides is governed by a trade-off between sensing and uptake"

See Methods for details.

This code was run on MATLAB R2020a, on a Mac computer.

Demo instructions:

Open and run the following .m files to obtain the results presented in the manuscript. Sample expected output is contained in the subfolders.
 
 - A_system_identification
	determines the parameter values that optimize the specified goodness of fit measure, using the model specified in line 155. (Functions of the different available models are in folder free_energy_models.

- B_SPECS_simulation
	runs an agent-based simulation of the cells using new transport-and-sensing model. It extends the SPECS simulation from Jiang, et al. and was used to ensure that the analytical approximation used in the system identification well matched the mechanistic model used in the simulations.


Available upon request:
  - Raw video of motile bacteria
  - Raw bacterial position data 
  - MATLAB code used to generate all plots



