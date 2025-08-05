/*------------------------------------------------------------------------*/
/*! \file poly_parser.h
    \brief contains helper functions for reading polynomials and substitution
    steps from file input.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

#ifndef POLY_PARSER_H_
#define POLY_PARSER_H_

// std includes.
#include <stdlib.h>

#include <fstream>
#include <regex>

// Local includes.
#include "polynom.h"

/**
    Inititate gates and specification polynomial.

    @param filename Name of the file with specification descriptions.

    @return Polynomial for specification
*/
void init_spec(Polynom & spec, std::string filename);

/**
    Creates specification polynomial from string line.

    @param line input string with polynomial description

    @return Polynomial for specification
*/
void read_spec_poly(Polynom & spec, std::string line);

/**
    Reduce spec by polynomials given in the file.

    @param spec specification Polynomial which will be reduced at the end
    
    @param gates array of gates which are used to save vars
    
    @param filename name of file containing the reduction polynomials

    @return reduced Polynomial
*/
void reduce_poly(Polynom & spec, std::string filename);

/**
    Convert given string line to a polynomial.
    
    @param spec specification Polynomial to be reduced

    @param line line description of a polynomial
*/
void reduce_by_one_line(Polynom & spec, std::string line); 

#endif /* POLY_PARSER_H_ */
