/*------------------------------------------------------------------------*/
/*! \file demo.cpp
    \brief demo file providing API examples for our package FastPoly.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/
#define VERSION "1.0"
/*------------------------------------------------------------------------*/

#include "poly_parser.h"
/*------------------------------------------------------------------------*/

/**
    Main Function of FastPoly Demo for different API usage.
    
*/
int main(int argc, char ** argv) {
  
  // API Example 1) Reading in from external file. 
  
  Polynom spec1;
  init_spec(spec1, "./demo/fulladder_example.txt");
  reduce_poly(spec1, "./demo/fulladder_example.txt");
  std::cout << "Result1:" << spec1 << std::endl;
  
  // API Example 2) Building polynomials from scratch.
  
  Polynom x8poly(6);  // Initialize poly for x8 gate. The initialization parameter has to be the maximum variable index which will be used in the polynomial.
  x8poly.createMonom(5); // Create and add monom 1*x5.
  x8poly.createMonom(6); // Create and add monom 1*x6.
  x8poly.createMonom(5,6, mpz_class(-1));  // Create and add monom -1*x5*x6.
  Polynom x7poly(4);  // Initialize poly for x7 gate.
  x7poly.createMonom(3); // Create and add monom 1*x5.
  x7poly.createMonom(4); // Create and add monom 1*x6.
  x7poly.createMonom(3,4, mpz_class(-2));   // Create and add monom -2*x3*x4.
  Polynom x6poly(4);  // Initialize poly for x6 gate.
  x6poly.createMonom(3,4);  // Create and add monom 1*x3*x4.
  Polynom x5poly(2);  // Initialize poly for x5 gate.
  x5poly.createMonom(1,2); // Create and add monom 1*x1*x2.
  Polynom x4poly(2);  // Initialize poly for x4 gate.
  x4poly.createMonom(1); // Create and add monom 1*x1.
  x4poly.createMonom(2); // Create and add monom 1*x2.
  x4poly.createMonom(1,2, mpz_class(-2));  // Create and add monom -2*x1*x2.
  Polynom spec2(8);  // Initialize spec2 polynomial.
  spec2.createMonom(8, mpz_class(2));  // Create monom 2*x8.
  spec2.createMonom(7);  // Create and add monom 1*x7.
  // Now replace x_i by corresponding polynomial x_ipoly.
  spec2.replaceVarByPoly(8, x8poly);
  spec2.replaceVarByPoly(7, x7poly);
  spec2.replaceVarByPoly(6, x6poly);
  spec2.replaceVarByPoly(5, x5poly);
  spec2.replaceVarByPoly(4, x4poly);
  std::cout << "Result2:" << spec2 << std::endl;
  
  // API Example 3) Using "shortcut" logic gate functions.
  Polynom spec3(8); // Initialize spec3 polynomial.
  spec3.createMonom(8, mpz_class(2));
  spec3.createMonom(7);
  spec3.replaceOR(8, 5, 6);
  // Function for replacing x8 by the OR gate polynomial with inputs x5 and x6.
  spec3.replaceXOR(7, 3, 4);
  spec3.replaceAND(6, 3, 4);
  spec3.replaceAND(5, 1, 2);
  spec3.replaceXOR(4, 1, 2);
  std::cout << "Result3:" << spec3 << std::endl;
  
  
  // Example 4) How to produce PAC proof files. 
  // NOTE: In the current version PAC proofs cannot be generated if phase optimization is applied inbetween reduction steps.
  Polynom spec4(8); // Inititate spec4 polynomial.
  spec4.createMonom(8, mpz_class(2));
  spec4.createMonom(7);
  // For starting proof generation, names of both output files have to be defined. 
  // The first name gives the output file where the starting polynomial as well as the polynomials
  // used for substitutions in every step are written to. The second name is the file for the actual PAC proof steps.
  spec4.startProofGeneration("example.polys", "example.proof");  
  spec4.replaceOR(8, 5, 6);
  spec4.replaceXOR(7, 3, 4);
  spec4.replaceAND(6, 3, 4);
  spec4.replaceAND(5, 1, 2);
  spec4.replaceXOR(4, 1, 2);
  std::cout << "Result4:" << spec4 << std::endl;
  
  // After all steps have been finished, the proof has to be generated with this function.
  // It is important to use the same file names as above when proof generation was started.
  writePolysIntoPACProof("example.polys", "example.proof");
  
  return 0;
}
