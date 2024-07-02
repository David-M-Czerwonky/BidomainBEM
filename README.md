Simulation code for Czerwonky, D. M., Aberra, A. S., & Gomez, L. J. (2023). A boundary element method of bidomain modeling for predicting cellular responses to electromagnetic fields. Journal of Neural Engineering. Read it here: https://iopscience.iop.org/article/10.1088/1741-2552/ad5704/meta  

Code Overview

We provide two folders of codes: 
  (1) the bidomain BEM codes
  (2) the hybrid cable codes

The BidomainCodes folder contains our implementations of bidomain integral equation derived in the above Journal of Neural Engineering Publication. This implementation leverages MATLAB with the support of C, and Fortran libaries to solve the bidomain integral equation. We provide testing scripts for simplistic transcanial electric stimulation (TES), transcranial magnetic stimulation (TMS), and deep brain stimulation (DBS) scenarios.  

The HybridCableCodes folder contains a MATLAB implementation of the Hybrid Cable approach as detailed in figure 1 of Joucla, S., Glière, A., & Yvert, B. (2014). Current approaches to model extracellular electrical neural microstimulation. Frontiers in computational neuroscience, 8, 13. The only difference is that instead of using an FEM solver, we use a modified set of 0th order BEM codes orginally from the following repository: https://github.com/luisgo/TMS_Efield_Solvers.


System requirements

All codes require a 64-bit system and Windows operating systems. All codes require an instillation of MATLAB 2021b or later. The bidomain solver codes are unacclerated so running the codes is computationally expensive. We recommend using simple test scenarios with less than 5,000 triangle elements for a system with 32 GB of RAM memory or less.
