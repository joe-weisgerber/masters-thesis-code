# masters-thesis-code

# Abstract

Driven by the growing sophistication of quantum simulators for studying lattice gauge theories, it has become essential to understand the dynamical properties of Abelian gauge theories using classical methods. One notable phenomenon is the occurrence of quantum many-body scars, which weakly violate the eigenstate thermalisation hypothesis. In this work, we employ tensor networks to study the dynamics of these states under various perturbations of the Hamiltonian. By linking the resulting thermalisation behaviour to whether a global, a local, or no symmetry is broken, we identify distinct time-evolution signatures that can be used to pinpoint specific perturbations in an experiment. Moreover, these dynamics also present an interesting theoretical question as they could deliver insights into the limitations of the eigenstate thermalisation hypothesis. While the full theoretical implications remain an active area of research, we provide different diagnostic methods for experimental setups and highlight the outstanding question of how spectral symmetry could provide a mechanism to drive highly oscillatory dynamics at early times.

# Code

The repository contains both the exact diagonalisation and the tensor network codes. The former are written in Python and the latter in Julia and make use of the ITensors.jl framework. 
Typically, we refer to the potential as p, charge conjugation breaking as g, constraint breaking as cb and gauge breaking as gb.
Some functionalities may not work as they have not been used at later stages of the project and have therefore not been kept up to date. 
We provide an example of how to use the exact diagonalisation code. The tensor network folder also contains the simulation outputs for L=16 as an h5 file, which may take some time to run otherwise.
