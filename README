Video presentation:
https://youtu.be/m5GvO2jOCF4

WARNING: This repository is based in the implementation of a NVIDIA GeForce 940MX GPU. If other GPU is employed, the architecture flag given to the compiler nvcc must be changed. (-arch=sm_50)

This repository contains code to implement the Boundary Wall Method in order to solve the 
Lippmann-Schwinger equation, which is an integral equation used in scattering theory to 
describe scattering processes given some potentials, which in this case, will be potentials
defined as contour integrals over parabollic contours weight by some gamma factor which
dependes on the parameter of the parabollical contour.

Analytical solutions are known for the problem of parabollic contours (see the jupyter notebook
LS.ipynb) and so this BWM can be proved to be efficient and approximate to analytical solutions.
By the way the Boundary Wall Method is based in matrices operations and many Hankel functions
evaluations in order to compute the scattered wave, parallel executions are proposed in order to
attack this problem and gain efficiency, so, CUDA and OMP parallel will be studied in this problem to determine some parallel metrics based on wall times and comparing this to the serial
implementation in C++.

The repository is organized as follows:

- Makefile: This contains every make command defined in order to execute different tasks which will be explained later.

- include/: Contains the main code files both in C++ and CUDA. For both of this it is
included the generation of parabolic billiards and the Boundary Wall Method.

- ComputationalTimes/: Contains the Python, C++ and CUDA files needed to implement the
time metrics and grafication of this for both C++ and CUDA implementations as matrix dimensions are bigger and optimizers for both compilers are included in the compilation flags.

- plot_dprob_dfase/: Contains the C++ and CUDA files for computing the scattered wave given the billiard parameters, the wave number and the angle of the incident wave in degrees. The plotting is done by matplotlib yn Python.

- plot_spectrum/: This directory contains the C++ code that calculates the resonances for this problem in a billiard. CUDA code is not included by the way this problem isn't really hard to compute.

- Videos/: This includes the C++ code that computes and organizes the different .dat outputs as the wave number or the angle of the inciden wave is changed. The animation is done in Python.

- LS.ipynb: This jupyter notebook contains the analytical solution for the parabollic problem (not confocal).


The makefile targets are the following:

- plotb: Plots the contour.

- plot: Plots the probability density and the phase of the scattered wave.

- plotspec: Plot the spectrum with the resonances identified.

- run_optimized: Executes the C++ and CUDA implementations with optimizers and registers the times recorded for different contour and grid number of points.

- run_not_optimized: Executes the C++ and CUDA implementations without optimizers and registers the times recorded for different contour and grid number of points.

- omp_times: Executes the C++ implementation for different number threads

- animation_k: Executes the multiple calls using parallel to compute the scattered wave by changing the wave number.

- animation_angles: Executes the multiple calls using parallel to compute the scattered wave as the angle of the incident wave is changed from 0 to 90 degrees.

WARNINGS: Before executing any task that saves data in .txt files be sure to execute "make clean" before doing so, in order to avoid accumulating data.