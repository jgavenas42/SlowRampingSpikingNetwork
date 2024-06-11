# SlowRampingSpikingNetwork
Code Repository for simulation and analysis of networks that exhibit slow-ramping activity as shown in Gavenas, Schurger, Rutishauser, & Maoz (2024) "Slow ramping emerges from spontaneous fluctuations in spiking neural networks."

snap_fluxsim.m runs the main simulations using the connMatrix.m function (for construction of connectivity matrix) and snap_singletrial.m function (for simulation of a single trial of spiking data). For simpler examples and to get a sense of how these parameters relate to network activity use snap_SimulateExample.m to simulate single trials and visualize activity.


Guide for plotting figures related to simulation data from Gavenas et al. (2024):

Figure 1 -- example activity in networks with different synaptic dynamics and clustering coefficients
1. Simulate and plot data using the MATLAB script snap_SimualteExample.m

Figure 2 & Figure S1 -- temporal characteristics of spontaneous fluctuations
1. Make sure you have the data fluxquant.csv and powspectrms.mat downloaded (should be in the fluxQuant folder).
	- Alternatively, you can simulate data using the MATLAB script snap_fluxSim.m and then calculate these statistics from simulated data by running the script snap_fluxcharacteristics.m
2. Code for generating figures is in the python notebook snap_fluctuationanalysis.ipynb and below the calculation of stats in the MATLAB script snap_fluxcharacteristics.m

Figure 3A, 3D left, 4C top, S2B, S2C -- comparison of slow ramping in networks with 100% fast and 50% slow synapses
1. Make sure you have the data downloaded. These data are in the folder data/flux
	- Alternatively, you can simulate data using the MATLAB script snap_fluxSim.m, just need to change some parameters including the Ksyns to 0, 0.5, and 1 from just 0 and 0.5.
2. Run the MATLAB script snap_PlotsMain.m 
3. The first sections of this script pertain to plotting the threshold-aligned results from simulated data.

Figure 3B, 4C bottom -- comparison of ramping in networks with varying synaptic ratios
1. Make sure you have the data downloaded. These data are in the folder data/Ksyn_cf
	- Alternatively, you can simulate data using the MATLAB script snap_KsynSim.m
2. Run the MATLAB script snap_ksynCF.m

Figure 5B -- pairwise correlations in neurons 
1. Make sure you have the data downloaded, they are in the folder pwCorr/flux/smoothedCorr_norm
2. Run the python notebook snap_pwCorrLME.ipynb

Figure S4D -- comparison of temporal generalization with ramping from "ground truth" slow ramping neurons.
1. To simulate and plot these data you should use the python notebook snap_poissonRamp.ipynb
