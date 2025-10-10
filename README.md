# Rt-cases-seqs
Rt from sequences vs case counts: simulation-based comparisons

This repo contains code and simulated outbreaks relating to the preprint *Comparing
methods to estimate time-varying reproduction numbers using genomic
and epidemiological data*, available at 
https://www.medrxiv.org/content/10.1101/2025.09.25.25336592v1.

## Code
In the codes folder, you'll find the following R files: 

auto_simulator.R: Generates outbreak data (linelist and FASTA files). You can specify
the number of outbreaks and other OOpidemic input parameters.

auto_practRt.R: Calculates the true Rt from the linelist.

auto_genTime_Incid.R: Calculates generation times and incidence from
the linelist.

auto_EpiEstimRt.R: Calculates the EpiEstim Rt.

auto_modfy_BST_ESTIM.R: Modifies BEAST’s Rt estimates.

auto_downSamp_fixed.R and auto_downSamp_switch.R: Used for
downsampling.

auto_process_plots.R: Post-processes results and generates plots.

building-intuition.R: Helps build intuition for understanding how DTW
works

In the main folder there are two R files, sequences-beast-setup.R and sequences-beast-collectRt.R, which set up
and then collect (respectively) the BEAST-based Rt estimates using the simulated sequences and the template.xml file. 

## Outputs
In the output folder there are folders with csv files corresponding to
simulated outbreaks in the paper. 

For example, EPiEstim_Rt contains 20 csv files corresponding to Rt distribution
estimates (ie with quantiles), through time, in simulated outbreaks.

Folders n_to_m_EPiEstim_Rt contain csv files corresponding to Rt distribution
estimates, through time, in simulated outbreaks, without or with
downsampling.

Other folders are similar. 

## Figures 

DTW_Alignment_Plots.pdf has examples to illustrate the DTW plot
approach (Figure 2 of the preprint). 

### Guide to figures in the paper and their filenames

Scenario A: optimum_combined_5plots

Scenario B: 70perc_combined_5plots

Scenario C: 10_combined_5plots

Scenario D: 10_to_70_combined_5plots

Scenario E: 70_to_10_combined_5plots

Figure 8: rsme_combined_plots (bar plots)




