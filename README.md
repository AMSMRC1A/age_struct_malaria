# age_struct_malaria

NB Some results are saved to named directories so the codebase folder should be downloaded in full and it's structure preserved so as not to introduce bugs.

Summary of key files:

Malaria_parameters_baseline.m (script) -> sets the basic model parameters as global variables. This file calls Malaria_parameters_transform.m to find the stable age 					      					      distribution via calls to balance_fertility.m to find a balanced birth rate function, if this option is selected, and then a call to find_stable_age.m

age_structured_Malaria.m (function) -> time-stepping algorithm for solving the ODE-PDE system. Calls Malaria_parameters_baseline.m to set parameters and find stable age distribution. Calls Malaria_IC and Immunity IC for initial/boundary conditions

run_age_structured_Malaria.m (script) -> time-stepping algorithm for solving the ODE-PDE system. Calls Malaria_parameters_baseline.m to set parameters and find stable age distribution. Calls Malaria_IC and Immunity IC for initial/boundary conditions

bifurcation_calcs.m (script) -> runs continuation algorithm tracking stability of the DFE and single endemic equilibrium versus chosen model parameter. Calls Malaria_parameters_baseline.m to set parameters and find stable age distribution

run_SA.m (script) -> runs sensitivity analysis of the model against each of a list of selected parameters. Calls Malaria_parameters_baseline.m to set the parameters [currently has some bugs]

run_age_structured_Malaria_parameterSearch.m (script) -> version of run_age_structured_Malaria.m to search the parameter space and assess impact of parameter changes in the sigmoids psi, phi and rho on the immunity distributions and endemic equilibrium














