# Rt-cases-seqs
Rt from sequences vs case counts: simulation-based comparisons

This repo contains code and notes relating to the preprint Comparing
methods to estimate time-varying reproduction numbers using genomic
and epidemiological data, available at 
https://www.medrxiv.org/content/10.1101/2025.09.25.25336592v1.

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

In the output folder there are folders with csv files corresponding to
simulated outbreaks in the paper. 

For example, EPiEstim_Rt contains 20 csv files corresponding to Rt distribution
estimates (ie with quantiles), through time, in simulated outbreaks.

Folders n_to_m_EPiEstim_Rt contain csv files corresponding to Rt distribution
estimates, through time, in simulated outbreaks, without or with
downsampling.

Other folders are similar. 

DTW_Alignment_Plots.pdf (and ..XXX.pdf, too) has examples to illustrate the DTW plot
approach (Figure 2 of the preprint). 




