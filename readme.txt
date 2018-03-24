%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%README FILE FOR "FIGURING OUT THE FED"
%BY CHRISTIAN MATTHES
%CONTACT: christian.matthes@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The main file is called main_clean.m. Running this file solves the model (using the posterior mode) and reproduces the main figures from the paper.
Auxiliary files include Paul Soderlind's code to solve for the optimal policy under discretion and commitment. 
The main file computes the error terms which can then be used to calculate the likelihood and carry out Bayesian inference.
The data is saved as both a .mat file and in ascii format (in that case every variable is saved in a separate file).