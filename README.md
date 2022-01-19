# RCT Simulator

A series of R files to simulate the results of multiple randomised trials, calculate Bayesian summary statistics, and display these graphically.

Contents:
* generator.R  
Generates datasets. The datasets used in the paper are in the /data/raw directory; you can generate your own (with different assumptions) by running generator.R - this is computationally demanding for large numbers of simulations.
* analysis.R  
Takes the datasets from the /data/raw directory and calculates Bayesian summary statistics and other summary values and exports the relevant dataframes, and saves them to the /data directory. These files are not in the github repo but are easily calculated from the raw data.
* output/outputs.rmd  
Takes the datasets produced by analysis.R and produces tables and figures.