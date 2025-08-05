/*------------------------------------------------------------------------*/
/*! \file proof_writer.h
    \brief contains helper functions for generating PAC proofs.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

#ifndef PROOF_WRITER_H_
#define PROOF_WRITER_H_

// std includes.
#include <stdlib.h>

#include <gmpxx.h>

#include <fstream>
#include <regex>
#include <string>

// Local includes.
//#include "polynom.h"
class Polynom;

//****************************************************************************************/
// Global variables
extern int axiomNum;
extern std::string polyfilename;
extern std::string prooffilename;
extern mpz_class modCoefProof;
extern bool firstline;


//****************************************************************************************/
/** Set the filenames where proof components will be written to.

	@param polyname std::string
	@param proofname std::string
*/	
void set_proof_filenames(std::string polyname, std::string proofname); 

/** Write starting polynomial as well as mod reduction coefficient and max. variable index of polynomial to proof file.

	@param inputpair std::pair<std::string, std::string> first is the polynomial string, second is the mod reduction number as string
	@param maxVarIndex int 
*/
void writeStartPolyToFile(std::pair<std::string, std::string> inputpair, int maxVarIndex);

/** Write the given string as PAC format axiom into polyfile.

	@param axiomStr std::string
*/
void writeNewPolyAxiom(std::string axiomStr);

/** Convert argument string to a polynomial axiom in PAC format.

	@param axiomStr std::string
*/
std::string convertPolyStringToPACFormat(std::string subStr);

/** Take a series of polynomial reduction steps given in inputName and write a complete PAC proof into outputName. 

	@param inputName std::string
	@param outputName std::string
*/
void writePolysIntoPACProof(std::string inputName, std::string outputName);

/** Initialize specification polynomial from file. 

	@param spec Polynom
	@param filename std::string
*/
void init_spec_from_PAC(Polynom& spec, std::string filename);

/** Read and create specification polynomial from string. 

	@param spec Polynom
	@param line std::string
*/
void read_spec_poly_from_PAC(Polynom& spec, std::string line);

/** Reduce specification polynomial by given series of polynomial reduction steps in file inputname
	and write PAC proof steps into file outputname.

	@param spec Polynom
	@param inputname std::string
	@param outputname std::string
*/
void reduce_poly_with_proof(Polynom& spec, std::string inputname, std::string outputname);

/** Reduce specification polynomial by one line representing a single polynomial reduction step.

	@param spec Polynom
	@param line std::string
	@param outputname std::string
	@param lineNum int
*/
void reduce_by_one_line_with_proof(Polynom& spec, std::string line, std::string outputname, int lineNum); 

/** Write the proof step of one step of polynomial reduction.

	@param outputname std::string
	@param usedAxiom int
	@param pol Polynom
	@param quotientStrVec std::vector<std::string>
*/
void writeOneLineIntoProof(std::string outputname, int usedAxiom, Polynom& pol, std::vector<std::string>& quotientStrVec);

/** Add a modulo reduction step to the proof step contained in currStr.

	@param pol Polynom
	@param currStr std::string
	@param modReduction mpz_class
*/
std::string addModReductionStep(Polynom& poly, std::string currStr, mpz_class modReduction);

/** Remove semicolons and leading line numbers from line.

	@param line std::string
*/
void removeLineNumAndSemicolon(std::string& line);


#endif /* PROOF_WRITER_H_ */
