ILQR Swarm Control using Random Finite Set Theory (Clohessy-Wiltshire Relative Motion Model) 
Version 1.00 10/2/2018

Algorithm used in the paper "Random Finite Set Thoery and Optimal Control for Large Spacecraft Swarms"
by: Bryce Doerr, Richard Linares, Pingping Zhu, Silvia Ferrari

MATLAB Code by Bryce Doerr doerr024(at)umn.edu with:
-GIF Plotting by Nicolae CINDEA from https://www.mathworks.com/matlabcentral/fileexchange/17463-movie-to-gif-converter?s_tid=prof_contriblnk
-Adapted code from Polygon point Plotting by Sulimon Sattari from https://ww2.mathworks.cn/matlabcentral/fileexchange/41454-grid-of-points-within-a-polygon

ILQR algorithm is based on Yuval Tassa's implementation
https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization
Tassa, Yuval, Tom Erez, and Emanuel Todorov. "Synthesis and stabilization of complex behaviors through online trajectory optimization." Intelligent Robots and Systems (IROS), 2012 IEEE/RSJ International Conference on. IEEE, 2012.

RFS Formulation based on Ba-Ngu Vo
Vo, Ba-Ngu, and Wing-Kin Ma. "The Gaussian mixture probability hypothesis density filter." IEEE Transactions on signal processing 54.11 (2006): 4091.

This is an implementation of Swarm ILQR control using Random Finite Set Theory.
The files are implemented as MATLAB scripts with the main file "main_ilqr.m"
The functions are commented for understanding.
This script has been successfully tested using the distance cost_nl_l2q.m. The gradients and hessians for cost_nl_l2q are in dj1.m file.
To run this script using cost_nl_cs or cost_nl_l2, the gradients and hessians must be determined. This is future work.
