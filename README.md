# MonteCarloHardDisk - updated April 6 2021 
This is a collaboration between the Cafiso lab and the DuBay lab at the University of Virginia. This project was published at https://pubs.acs.org/doi/abs/10.1021/jacs.0c01754

This model involved using Monte Carlo simulations on circular particles with two sticky interfaces. The interfaces were placed 167 degrees from each other, and we then analyzed how the radial distribution function changes with different interaction strengths and turning certain interfaces on. We collected RDFs for various spin sites (tags) that were indicated by the Cafiso group.

Folder Structure: 

/1_RunSimulation contains the scripts necessary to run a Monte Carlo simulation of the model. 

/2_CheckProgress contains the scripts used to moniter acceptance criteria 

/3_RDF is the analysis script to produce the radial distribution function. Here there is a subdirectory containing the text files from our simulations that are utilized within the paper 

/4_Snapshots scripts for visualization 
