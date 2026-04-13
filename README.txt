This zip file contains the data, simulation results, and R code needed to reproduce all tables and figures found in the paper titled "Experimental Designs for Electric Bus Routes."

Before running any of these files, you need to set the working directory to the root directory of this file (e.g. "C:/Users/riosn/OneDrive/Desktop/DASHBusRCodeR1"). This can be done in RStudio by going to Session -> Set Working Directory and navigating to the correct folder. This only needs to be done ONCE - all other files use relative paths to this working directory. If a file changes the working directory to another directory, it navigates back to the root directory afterwards (to save you time). 

Here is a description of all folders and files:

- SimResults: This folder contains simulation results in csv files that are needed to reproduce Tables 1 and 2 from Section 3 of the paper.

- bestSAdesign_n50L5.csv, bestSAdesign_n60L5.csv: These files store the best designs found in Appendix C, Tables C1 and C2, respectively.

- bus_data_X_Y_format.rds: This is an RDS object that stores the DASH bus data that are used throughout Section 4.

- bus_graph_rds: This is an RDS object that stores the bus network for the DASH bus data.

- Edge_prior.RData: This RData object stores information used to find the distance-based prior (mu2) in Section 4. This file is required to run Section4.R.

- Figures1and2.R: Running this R code produces Figure 1 and Figure 2.

- Figure3.R: Running this R code produces Figure 3.

- SA_rel_Deffs_n50_L5.csv, SA_rel_Deffs_n50_L6.csv: These csv files can also be obtained by running Section5.R. They are called by Figure3.R to produce Figure 3.

- Section4.1Sim.R: This R code was used to run the numerical studies in Section 4.1.

- Section4.2Sim.R: This R code was used to run the numerical studies in Section 4.2.

- Section5.R: This R code reproduces Table 3, Table 4, and Tables C1 and C2. It essentially runs all of the analysis done in Section 5.

- source_code.R: This R code contains custom functions that are used throughout the R files here. 

- Table1.R: Running this R code reproduces Table 1.

- Table2.R: Running this R code reproduces Table 2.

- TableB1.R: Running this R code reproduces Table B1 in Appendix B. 

- TableB2.R: Running this R code reproduces Table B2 in Appendix B.

- TableB3.R: Running this R code reproduces Table B3 in Appendix B.

- AppendixBDesignSim.R: This code was used to run the simulations to make Tables B1 and B2.

- AppendixBSimTableB3.R: This code was used to run the simulations used to make Table B3.

- AppendixD.R: This code reproduces the analysis in Appendix D. It will reproduce Table D1 and D2.