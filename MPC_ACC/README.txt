MPC Swarm Control using Random Finite Set Theory (Acceleration Model)
Version 1.00 10/2/2018

Algorithm used in the paper "Random Finite Set Theory and Optimal Control for Large Spacecraft Swarms"
by: Bryce Doerr, Richard Linares, Pingping Zhu, Silvia Ferrari

MATLAB Code by Bryce Doerr doerr024(at)umn.edu with:
-uses MATLAB's fminunc, therefore you need the Optimization Toolbox

RFS Formulation based on Ba-Ngu Vo
Vo, Ba-Ngu, and Wing-Kin Ma. "The Gaussian mixture probability hypothesis density filter." IEEE Transactions on signal processing 54.11 (2006): 4091.

This is an implementation of Swarm MPC control using Random Finite Set Theory.
The MPC control uses the fminunc solver and the Quasi-Newton algorithm.

The files are implemented as MATLAB scripts with the main file "main_mpc.m"
The functions are commented for understanding.
This script has been successfully tested using all three distances in the paper. The CS divergence does not converge for this problem.
