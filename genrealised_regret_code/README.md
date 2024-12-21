# Generalised Regret Optimal Controller Synthesis for Constrained Systems
Code and simulation results accompanying the paper "Generalised Regret Optimal Controller Synthesis for Constrained Systems" by Alexandre Didier ([ORCID](https://orcid.org/0000-0001-5396-8996)) and Melanie N. Zeilinger ([ORCID](https://orcid.org/0000-0003-4570-7571)), published at the IFAC World Congress 2023. 
Created by Alexandre Didier, Prof. Dr. Melanie N. Zeilinger, Institute for Dynamic Systems and Control, ETH Zurich. March 2023.

### Files
main.m: main executable Matlab file\
cleanup.m: file for cleaning up after code is run\
Controller.m: abstract controller class\
setup.m: file for setting up the environment\
sim_data.mat: saved simulation data, can be loaded with Matlab

controllers/GenReg_constr.m: constrained generalised regret controller class\
controllers/GenReg_PWB_constr.m: constrained generalised regret controller class with pointwise bounded disturbances\
controllers/H_2_constr.m: constrained H_2 controller class\
controllers/H_inf_constr.m: constrained H_infinity controller class\
controllers/NC_constr.m: constrained non causal controller class

noise/constant.m: constant noise function\
noise/gaussian_trunc_ell.m: truncated Gaussian noise on an ellipsoid function\
noise/sawt.m: sawtooth noise function\
noise/sine.m: sine noise function\
noise/stair.m: stair noise function\
noise/step.m: step noise function\
noise/uniform_ell.m: uniform noise function on an ellipsoid

params/get_params.m: loads simulation parameters\
params/params_rocket.m: containts simulation parameters

sim/sim_rocket.m: simulates rocket with considered controllers and noise profiles

### License
BSD-2 Clause License. See license file for more information

### Getting Started
Install [Matlab](https://ch.mathworks.com/products/matlab.html) (used version: R2020b), [YALMIP](https://yalmip.github.io/) (used version: 20210331) and possibly [MOSEK](https://www.mosek.com/) (used version: 9.2.29) to solve the provided SDPs (you can change the SDP solver in sim/sim_rocket.m).\
Machine used has an Intel(R) Core(TM) i7-10510U CPU.

### Contact
To contact me, please send an email to adidier@ethz.ch.
