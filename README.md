# Bruggeman_Effective_Medium

Description:                                                                                                                                               
============
 
Bruggeman_Effective_Medium is a collection of C++ codes
for the computation of an effective complex index of
refraction based on the approximate formulae published
by Bruggeman (1935).


Getting Started:
================

The first step consists of the computation of the Bruggeman
effective refractive index using ./get_index_exec.
The second step is the computation of a real part of the
refractive index that fulfills the Kramers-Kronig relation
by executing ./KK_analysis_exe.

Requirements:
=============

The code requires boost, version 1.63.0 or higher and the GNU Scientific Library.


Compiling the code:
===================

The code for the effective medium calculations is compiled using the `makefile` file, while the subsequent code for the Kramers-Kronig integral is compiled using `makefile_KK`.
The compilation has been tested with LLVM clang++ and gcc version 9.1.0.
You may need to update the compiler and the library include directories in the provided makefiles prior to a successful compilation.


Using the Code:
===============

- If you are using the code, please cite the corresponding paper by Stegmann and Yang, (2017) in any publication (journal paper, presentation, poster,...) where you are using the code's results.
- If you plan to develop or modify the code, please create a `feature` branch from `develop` or fork the repository.


LICENSE:
========

Please see the file LICENSE.txt for details
 

References:
===========

**P. G. Stegmann, and P. Yang**, (2017). *A regional, size-dependent, and causal effective medium model for Asian and Saharan mineral dust refractive index spectra*, J. Aer. Sci. 114, pp. 327-341.
 
**Bruggeman, D. A. G.**, (1935). *Berechnung verschiedener physikalischer Konstanten von heterogenen Substanzen: 1. Dielektrizitätskonstanten und Wärmeleitfähigkeiten der Mischkörper aus isotropen Substanzen*, Annalen der Physik, 5. Folge, Band 24, pp. 636-664. Transl.: Calculation of various physical constants of     heterogeneous substances: 1. Dielectric permittivities and heat conductivities of mixed bodies made out of isotropic substances.
 
