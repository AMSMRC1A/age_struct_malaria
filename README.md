# Project Team

This repository is comprised of codes related to the project "Modeling Immunity to Malaria with an Age-Structured PDE Framework" arising from the American Mathematical Society's Mathematical Research Community on infectious diseases (2020-2021). The team members on this project are:

Lauren Childs (Virginia Tech), Christina Edholm (Scripps College), Denis Patterson (Princeton University), Joan Ponce (UCLA), Olivia Propser (University of Tennessee, Knoxville), Zhuolin Qu (University of Texas at San Antonio) and Lihong Zhao (University of California Merced)

# ReadMe

I want to:

1. Run a single simulation -> run_age_structured_Malaria.m 
2. Run bifurcation analysis -> run_bifurcation_calcs.m
3. Run sensitivity analysis -> run_SA.m
4. Calibrate sigmoids -> run_parameter_search.m -> Data_Fitting/run_parameter_fit.m

NB Some results are saved to named directories so the codebase folder should be downloaded in full and it's structure preserved so as not to introduce bugs.

# Summary of key files

Malaria_parameters_baseline.m (script) -> sets the basic model parameters as global variables. This file calls Malaria_parameters_transform.m to find the stable age 					      					      distribution via calls to balance_fertility.m to find a balanced birth rate function, if this option is selected, and then a call to find_stable_age.m

age_structured_Malaria.m (function) -> time-stepping algorithm for solving the ODE-PDE system

run_age_structured_Malaria.m (script) -> Calls age_structured_Malaria.m for time-stepping algorithm to solve the ODE-PDE system. Calls Malaria_parameters_baseline.m to set parameters and find stable age distribution. Calls Malaria_IC and Immunity IC for initial/boundary conditions

bifurcation_calcs.m (script) -> runs continuation algorithm tracking stability of the DFE and single endemic equilibrium versus chosen model parameter. Calls Malaria_parameters_baseline.m to set parameters and find stable age distribution

run_SA.m (script) -> runs sensitivity analysis of the model against each of a list of selected parameters. Calls Malaria_parameters_baseline.m to set the parameters [currently has some bugs]

run_parameter_search.m (script) -> version of run_age_structured_Malaria.m to search the parameter space and assess impact of parameter changes in the sigmoids psi, phi and rho on the immunity distributions and endemic equilibrium

Data_Fitting (folder) -> Data sets from Filipe and R-B papers for fitting the sigmoids. 














