# OCDES: OCDE Simulator
Xiao Zhao, Forschungszentrum Jülich, 07.2020

This code is an extension of my work at Forschungszentrum Jülich. 

# Overview
OCDES is a MATLAB-based tool that performs numerical integration to solve Optimization-Constrained Differential Equations (OCDE):                                                            
dx = f(x,v), x(0)=x_0,	 (1a)

v∈arg min_v ⁡g(x,v),	    (1b)

s.t.h_i(x,v)=0, i=1,…,M, (1c)

   l_j(x,v)≥0,  j=1,…,N. (1d)

x∈R^m and v∈R^n. f:R^m×R^n→R^m, g:R^m×R^n→R,h_i:R^m×R^n→R , l_j:R^m×R^n→R are at least second order continuously differentiable. OCDES provides efficient numerical solution of OCDE by using local optimality condition. A sequence of DAE systems are generated and classical index-1 DAE simulator is applied to solve the derived DAE systems.

# Basic Requirements

Matlab, version 2014 or higher.

Matlab Symbolic Toolbox

# Installation
The simulator needs Matlab environment. 

# How to use OCDES
Please refer to the example demo1.m. General steps of using OCDES are:

(1) Define state variables x and optimization variables v in symbolic format. 

(2) Define functions f,  g, h, l in symbolic format.

(3) Give initial condition x(0) and initial guess of v(0).

(4) Specify options for solving the inner NLP for initialization

opt_init.tol_act:	tolerance to check active inequality constraints

opt_init.optimoptions:	Optimization options, cf. MATLAB optimoptions

opt_sol.MaxNoUptActiveSet:	maximum number of updating active set

opt_sol.tol_feasible:	feasibility tolerance

(5) Specify options for integration

tstart:	starting time of simulation

tfinal: 	ending time of simulation

opt_sol.integrator: Selected integration

opt_sol.opt_integrator:	Integration options 

(7) Call sOCDE_main.m to solve the OCDE.

# Citation
Please cite [1], if you use the code.

# License
Copyright (c) 2020: Forschungszentrum Jülich GmbH, Jülich, Germany. 

Author: Xiao Zhao, Email: Xiao.Zhao@outlook.de

# References
[1] Zhao, Ploch, Noack, Wiechert, Mitsos, von Lieres, Analysis of local well-posedness of optimization-constrained differential equations by local optimality conditions, AIChE J., DOI:10.1002/aic.16548.

