[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## FastPoly: An Efficient Polynomial Package for Symbolic Computer Algebra

by Alexander Konrad and Christoph Scholl, University of Freiburg

For further information we refer to the paper

Alexander Konrad, Christoph Scholl. 
 [`FastPoly: An Efficient Polynomial Package for Symbolic Computer Algebra.`]
In Proc. 25th Intl. Conf. on Formal Methods in Computer Aided Design (FMCAD'25), p. TBA, 2025.

  
----------------------------------------------------------------  
  
Dependencies: `libgmp` (https://gmplib.org/)

To use the provided polynomial package the GNU Multiprecision (GMP) library has to be installed first. 
It can either be aquired from https://gmplib.org/
or be installed directly on Ubuntu by the command:

sudo apt-get install libgmp-dev

----------------------------------------------------------------

Using the package:

To use the package in your own software, include the src folder in your source files 
(and compile them with your other source files)
and include the poly_parser.h in your code (like it is done in the demo).

----------------------------------------------------------------

Provided demo:

We provide a demo file fastpoly_demo.cpp which shows how to use the different APIs of
our package. Use `./configure.sh && make` to configure and build the demo. 

To run the demo: `./fastpoly_demo`

Example 1) 	Shows how to initialize a polynomial from a file and perform substitution steps written in the same file.

Example 2) 	Shows how to initialize a starting polynomial as well as gate polynomials using the constructors.
		Substitution steps are invoked by the "replaceVarbyPoly()" function.

Example 3)	Shows how to use "shortcut" functions for replacing variables by functions of predefined basic logic gates.

Example 4)	Shows how to start proof generation.

----------------------------------------------------------------

The Benchmark folder provides benchmarks used in the paper stated above. 
For details on the benchmarks we refer to the paper.

----------------------------------------------------------------
