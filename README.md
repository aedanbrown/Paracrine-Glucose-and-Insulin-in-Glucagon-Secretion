# Paracrine Glucose and Insulin in Glucagon Secretion
 This repository contains the code used to develop the figures in the paper.

## Running the code
 The code used to generate the data in Figures 3 through 5 and Supplemental Figures 4 and 5 are in the Simulations folders. To run the code within these folders, code from the Base Models and Simulate Functions must be copied into the same folder as the main code. The table below shows the files that must be copied for each file in Curvefitting and Simulations.

|File| Description| Additional Files|
|--- | --- | --- |
|run_alphaBetaModel_batch_multg.m| Calculate total insulin and glucagon secretion in a batch setting - used for Figure 3 and Supplemental Figures 4. |simulate_alphaBetaModel_batch.m, alphaBetaModel_batch.m|
|run_alphaBetaModel_perfusion_multg.m| Calculate total insulin and glucagon secretion in a perfusion setting - used for Figures 3 and 5 and Supplemental Figures 4 and 5. | simulate_alphaBetaModel_perfusion.m, alphaBetaModel_perfusion.m|
|run_alphaBetaModel_perfusion_multg_invivo.m| Calculate total insulin and glucagon secretion in a in vivo perfusion setting for Figure 4. | simulate_alphaBetaModel_perfusion.m, alphaBetaModel_perfusion.m|
|perfusion_sensitivityAnalysis.m| Perform the sensitivity analysis to create the data for Supplemental Figure 3.| simulate_alphaBetaModel_perfusion.m, alphaBetaModel_perfusion.m|


**References cited in the code are in References.txt**
