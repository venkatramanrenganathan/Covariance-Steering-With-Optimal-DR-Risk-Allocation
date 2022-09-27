# Covariance Steering With Optimal Distributionally Robust (DR) Risk Allocation 
This repository contains the MATLAB code for simulating covariance steering for stochastic linear systems using optimal and distributionally robust risk allocation.

![covariance_steering](https://github.com/venkatramanrenganathan/Covariance-Steering-With-Optimal-DR-Risk-Allocation/blob/main/Figures/polyspaceDRIRA2Dtraj.png)

# Dependencies
- Matlab
- Yalmip
- MOSEK solver

# Procedure to run the code
1. To just run covariance steering with polytopic distributionally robust state risk constraints, run the matlab code `mainDRCS.m` which will load the required system data and generate the desired plots.
2. To run covariance steering with iterative risk allocation for polytopic distributionally robust state risk constraints, run the matlab code `mainCSDRRA.m` which will load the required system data and generate the desired plots.
3. For conic state risk constraints, run the matlab code `mainCSDR_CC.m` which will load the required system data and generate the desired plots.

## Variations while running `mainDRCS.m`, `mainCSDRRA.m` and `mainDRCS_CC.m` files
    * Set `riskSelectFlag = 1` in line `13` for running simulations with Gaussian chance constraints
    * Set `riskSelectFlag = 2` in line `13` for running simulations with Distributionally Robust risk constraints 
    * Set `dynamicsSelectFlag = 1` in line `14` for running simulations with 3D spacecraft dynamics
    * Set `dynamicsSelectFlag = 2` in line `14` for running simulations with Double Integrator dynamics

# Contributing Authors
1. [Venkatraman Renganathan - Lund University](https://github.com/venkatramanrenganathan)
2. [Joshua Pilipovsky - Georgia Institute of Technology](https://github.com/JoshPilipovsky)
3. [Panagiotis Tsoitras - Georgia Institute of Technology](https://dcsl.gatech.edu/tsiotras.html)

# Funding Acknowledgement
1. V. Renganathan's work has been supported by the *European Research Council (ERC)* under the European Union’s Horizon 2020 research and innovation program under grant agreement No `834142 - Scalable Control`.
2. The work of the J. Pilipovsky and P. Tsoitras has been supported by `NASA University Leadership Initiative award 80NSSC20M0163` and `ONR award N00014-18-1-2828`.

# Affiliation
1. [Department of Automatic Control, Lund University, Sweden](https://control.lth.se)
2. [Dynamics and Control Systems Laboratory (DCSL) at Georgia Tech, USA](https://dcsl.gatech.edu)
