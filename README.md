# Paracrine_Interactions_In_Vitro
 This repository contains the code used to develop the figures in the paper.

## Running the code
 The code used to generate the data in Figures 2-6 and Supplemental Figures 1-3 is contained in the Curvefitting and Simulations folders. To run the code within these folders, code from the Base Models and Simulate Functions must be copied into the same folder as the main code. The table below shows the files that must be copied for each file in Curvefitting and Simulations.

|File| Description| Additional Files|
|--- | --- | --- |
|run_alphaBetaModel_batch_multg.m| Calculate total insulin and glucagon secretion in a batch setting - used for Figure 3A and Supplemental Figure 3A. |simulate_alphaBetaModel_batch.m, alphaBetaModel_batch.m|
|run_alphaBetaModel_perfusion_multg.m| Calculate total insulin and glucagon secretion in a perfusion setting - used for Figures 3, 4, and 5 and Supplemental Figure 3A. | simulate_alphaBetaModel_perfusion.m, alphaBetaModel_perfusion.m|
|run_alphaBetaModel_perfusion_multg_F6.m| Calculate total insulin and glucagon secretion in a perfusion setting using the parameters for Figure 6. | simulate_alphaBetaModel_perfusion.m, alphaBetaModel_perfusion.m|
|perfusion_sensitivityAnalysis.m| Perform the sensitivity analysis to create the data for Supplemental Figure 2| simulate_alphaBetaModel_perfusion.m, alphaBetaModel_perfusion.m|
|CF_steadyStateInsulinSecretion_human.m| Curve fit the steady-state insulin secretion in humans - Supplemental Figure 1B. | |
|CF_steadyStateInsulinSecretion_mice.m| Curve fit the steady-state insulin secretion in mice - Supplemental Figure 1A.| |
|compareCF_alphaBetaModel_batch_humans.m| Compare the results of the curve fit of human islets in a batch setting - Figure 2D. | simulate_alphaBetaModel_batch.m, alphaBetaModel_batch.m|
|compareCF_alphaBetaModel_batch_mice.m | Compare the results of the curve fit of mice islets in a batch setting - Figure 2B. | simulate_alphaBetaModel_batch.m, alphaBetaModel_batch.m|
|compareCF_kineticInsulinModel_perfusion_human.m| Compare the results of the curve fit of human islets in a perfusion setting - Figure 2C. | simulate_kineticInsulinModel_perfusion_multi.m, simulate_kineticInsulinModel_perfusion.m, kineticInsulinModel_perfusion.m|
|compareCF_kineticInsulinModel_perfusion_mice.m| Compare the results of the curve fit of mice islets in a perfusion setting - Figure 2A. | simulate_kineticInsulinModel_perfusion.m, kineticInsulinModel_perfusion.m|





**References cited in the code are in References.txt**