/*------------------------------------------------------------------------*/
/*! \file polynom.h
    \brief contains the class Polynom for representing polynomials and 
    defining all polynomial manipulation functions.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

#ifndef POLYNOM_H_
#define POLYNOM_H_

// std includes.
#include <stdlib.h>
#include <list>
#include <set>
#include <iostream>
#include <climits>
#include <vector>
#include <cassert>
#include <deque>
#include <cstdint>
#include <regex>

// Local includes.
#include "monom.h"
#include "proof_writer.h"

class Polynom {
	
	friend class Circuit;
	public:
		//*********************** Constructors  ******************************************************//
		
		/** Default Constructor. */
		Polynom();
		
		// Construct polynomial. varSize is the max. variable index which will be inserted.
		/** Constructor for polynomial. varSize is the max. variable index which will be inserted. 

			@param varSize int
		*/
		Polynom(int varSize); 
		
		/** Copy constructor.

			@param old Polynom to copy
		*/
		Polynom(const Polynom& old);
		
		/** Assignment operator.

			@param old Polynom to copy
		*/
		Polynom& operator=(const Polynom& old);
		
		/** Destructor. */
		virtual ~Polynom();
		
		//***************** Functions for finding, adding and removing parts of polynomial, mostly monomials  *********************//
		
		/** Add monomial to the polynomial.

			@param mon Monom
			@return Monom* pointer to just added monomial
		*/
		Monom* addMonom(Monom mon);
		
		/** Erase monomial from the polynomial.

			@param mon Monom
		*/
		void eraseMonom(Monom mon);
		
		/** Add reference variable to corresponding refList and monomials' pointers. 

			@param mon Monom
			@param index varIndex
			@param i int pos of variable in the monomial
		*/
		void addRefVar(Monom& mon, varIndex index, int i);
		
		/** Return vector of pointers to all monomials which contain the searched monomial.

			@param mon Monom searched for
			@return std::vector<Monom*> vector of pointers to all monomials containg mon.
		*/
		std::vector<Monom*> findContaining(Monom& mon);
		
		/** Return vector of pointers to all monomials which contain the searched variable.

			@param var varIndex searched for
			@return std::vector<Monom*> vector of pointers to all monomials containg var.
		*/
		std::vector<Monom*> findContainingVar(varIndex var);
		
		/** If mon contained in polynomial, return pointer to it. Otherwise return NULL.

			@param mon Monom searched for
			@return Monom*
		*/
		Monom* findExact(Monom& mon);
		
		/** Return true, if variable var is contained in polynomial.

			@param var varIndex searched for
			@return bool
		*/
		bool containsVar(varIndex var);
		
		/** Add complete polynomial to this polynomial. Return true if successful.
			If varSize of this polynomial is smaller than varSize of other, the polynomial
			cannot be added, returning false. 

			@param other Polynom
			@return bool
		*/
		bool addPolynom(const Polynom& other);
		
		/** Parse polynomial from string, adding all monomials to it.
			The polynomial already has to be initialized.

			@param std::string
		*/
		void parsePolyFromString(std::string);

		//******************************* API functions for directly creating monomials  ***********************************//
		
		/** Create and add monomial with size 1 and variable index1. Coefficient can be set directly.

			@param index1 varIndex
			@param coef mpz_class
		*/		
		void createMonom(varIndex index1, mpz_class coef=0);
		
		/** Create and add monomial with size 2 and variables index1 and index2. Coefficient can be set directly.

			@param index1 varIndex
			@param index2 varIndex
			@param coef mpz_class
		*/	
		void createMonom(varIndex index1, varIndex index2, mpz_class coef=0);
		
		/** Create and add monomial from integer array. Coefficient can be set directly.

			@param myints[] varIndex
			@param size int
			@param coef mpz_class
		*/
		void createMonom(varIndex myints[], int size, mpz_class coef=0);
		
		//******************************* Functions for controlling modulo reduction  ***********************************//
		
		/** Change mode of modulo reduction.

			@param mode bool
		*/
		void setModReduction(bool mode);
		
		/** Set modulo reduction coefficient. Every monomial coefficient will be reduced by it.

			@param modNum mpz_class
		*/
		void setModReductionNumber(mpz_class modNum);
		
		/** Reduce all monomial coefficients by modNum.

			@param modNum mpz_class
		*/
		void modReducePoly(mpz_class modNum);
		
		/** Reduce all monomial coefficients by modNum. Return the quotient from the corresponding polynomial division.

			@param modNum mpz_class
			@return std::vector<Monom> quotient of the polynomial division with modNum
		*/
		std::vector<Monom> modReductionWithQuotient(mpz_class modNum);
		
		/** Reduce all monomial coefficients by modNum. Return the quotient from the corresponding polynomial division as std::string.

			@param modNum mpz_class
			@return std::vector<std::string> quotient of the polynomial division with modNum as std::strings.
		*/
		std::vector<std::string> modReductionWithQuotientStr(mpz_class modNum);
		
		//******************************* Functions for manipulating phases of variables  ***********************************//
		
		/** Change phase of a variable by replacing it with its own negation.

			@param replace varIndex variable which phase will be flipped
		*/
		void negateVar(varIndex replace);
		
		/** Change phase of a variable. Improved implementation which uses some internal properties to improve the phase flipping.

			@param replace varIndex variable which phase will be flipped
		*/
		void negateVarImproved(varIndex var);
		
		/** Change phase of a variable. Improved implementation which uses some internal properties to improve the phase flipping.
			Additionally save the quotient of the "polynomial division" carried out by this negation.

			@param replace varIndex variable which phase will be flipped
			@param quotient std::vector<Monom> actually not used in current implementation
			@param quotientStrVec std::vector<std::string>
		*/
		void negateVarImprovedWithQuotient(varIndex var, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec);
		
		/** Return how much the polynomial size will change if the variable phase gets flipped in the given monomial. 
			This function only predicts how much the polynomial size will change, it does not apply the actual phase flipping.  

			@param mon Monom
			@param var varIndex
		*/
		int phaseChangeEffectOnMonom(Monom& mon, varIndex var);
		
		/** Apply a greedy phase change by flipping phases of all variables one after another, in ascending order. 
			If phase flipping a variable does not reduce polynomial size, it is reverted before continuing with next variables.   

			@return int by how much the polynomial size has improved
		*/
		int greedyPhaseChange();
		
		/** Again a greedy phase change, but in descending order.    

			@return int by how much the polynomial size has improved
		*/
		int greedyPhaseChangeBackward();
		
		/** Greedy phase change only of the variables given in the argument vector.    

			@param signalsToChange std::vector<varIndex>
			@return int by how much the polynomial size has improved
		*/
		int greedyPhaseChangeCustom(std::vector<varIndex>& signalsToChange);
		
		/** Greedy phase change only of the variables given in the argument list.    

			@param signalsToChange std::list<varIndex>
			@return int by how much the polynomial size has improved
		*/
		int greedyPhaseChangeCustom(std::list<varIndex>& signalsToChange);
		
		/** Greedy phase change only of the variables given in the first argument list.
			Additionally write changed phases into second argument list.   

			@param signalsToChange std::list<varIndex>
			@param changedPhases std::list<uint32_t>
			@return int by how much the polynomial size has improved
		*/
		int greedyPhaseChangeCustom(std::list<varIndex>& signalsToChange, std::list<uint32_t>& changedPhases);
		
		/** Change phase of variable. If polynomial size reduces, return true. If not, revert the phase flip and return false.     

			@param var varIndex
			@return bool
		*/
		bool testPhaseChangeSingleVariable(varIndex var);
		
		/** Improved phase changing of variable. If polynomial size reduces, return true. If not, revert the phase flip and return false.     

			@param var varIndex
			@return bool
		*/
		bool testPhaseChangeSingleVariableImproved(varIndex var);
		
		/** Print out all variables which phase is currently negated. */
		void reportVarPhases();
		
		//************** Functions for replacing one variable in the polynomial by gate functions.  *********************//

		/** Replace variable "replace" by the polynomial for a logic AND gate with inputs in1 and in2.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
		*/	
		void replaceAND(varIndex replace, varIndex in1, varIndex in2);
		
		/** Replace variable "replace" by the polynomial for a logic OR gate with inputs in1 and in2.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
		*/
		void replaceOR(varIndex replace, varIndex in1, varIndex in2);
		
		/** Replace variable "replace" by the polynomial for a logic XOR gate with inputs in1 and in2.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
		*/
		void replaceXOR(varIndex replace, varIndex in1, varIndex in2);
		
		/** Replace variable "replace" by the polynomial for a logic AND gate with inputs in1 and in2 where in1 is negated.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
		*/
		void replaceANDOneNegation(varIndex replace, varIndex in1, varIndex in2);
		
		/** Replace variable "replace" by the polynomial for a logic OR gate with inputs in1 and in2 where in1 is negated.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
		*/
		void replaceOROneNegation(varIndex replace, varIndex in1, varIndex in2);
		
		/** Replace variable "replace" by the polynomial for a logic XOR gate with inputs in1 and in2 where in1 is negated.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
		*/
		void replaceXOROneNegation(varIndex replace, varIndex in1, varIndex in2);
		
		/** Replace variable "replace" by the polynomial for a logic AND gate with inputs in1 and in2 where both inputs are negated.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
		*/
		void replaceANDDoubleNegation(varIndex replace, varIndex in1, varIndex in2);
		
		/** Replace variable "replace" by the polynomial for a logic OR gate with inputs in1 and in2 where both inputs are negated.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
		*/
		void replaceORDoubleNegation(varIndex replace, varIndex in1, varIndex in2);
		
		/** Replace variable "replace" by the polynomial for a logic NOT gate with input in1.  

			@param replace varIndex variable to replace
			@param in1 varIndex
		*/
		void replaceNOT(varIndex replace, varIndex in1);
		
		/** Replace variable "replace" by the polynomial for a logic BUFFER gate with input in1.  

			@param replace varIndex variable to replace
			@param in1 varIndex
		*/
		void replaceBUFFER(varIndex replace, varIndex in1);
		
		/** Replace variable "replace" by constant 0.  

			@param replace varIndex variable to replace
		*/
		void replaceCON0(varIndex replace);
		
		/** Replace variable "replace" by constant 1.  

			@param replace varIndex variable to replace
		*/
		void replaceCON1(varIndex replace);
		
		/** Replace variable "replace" by the polynomial for a logic AND gate with inputs in1 and in2. Negations on the inputs can
			be indicated by phase1 and phase2.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
			@param phase1 bool false indicates negation on in1
			@param phase2 bool false indicates negation on in2
		*/
		void replaceANDDependingOnNegations(varIndex replace, varIndex in1, varIndex in2, bool phase1, bool phase2);
		
		//***************** General function for replacing a variable in the polynomial.  *************************************//
		
		/** Replace variable "replace" by a list of monomials.  

			@param replace varIndex variable to replace
			@param mons std::list<Monom>
		*/
		void replaceVar(varIndex replace, std::list<Monom>& mons);
		
		/** Replace variable "replace" by a set of monomials.  

			@param replace varIndex variable to replace
			@param mons std::set<Monom>
		*/
		void replaceVar(varIndex replace, std::set<Monom>* mons);
		
		/** Replace variable "replace" by a polynomial.  

			@param replace varIndex variable to replace
			@param poly Polynom
		*/
		void replaceVarByPoly(varIndex replace, Polynom& poly);
		
		//************************ Replace functions with additional quotient return.  ***************************************//

		/** Replace variable "replace" by a list of monomials and additionally save the quotient of the "polynomial division" as vector of strings.

			@param replace varIndex variable to replace
			@param mons std::list<Monom>
			@param quotient std::vector<Monom> actually not used in current implementation
			@param quotientStrVec std::vector<std::string>
		*/
		void replaceVarWithQuotients(varIndex replace, std::list<Monom>& mons, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec);
		
		/** Replace variable "replace" by polynomial of logic AND gate and additionally save the quotient of the "polynomial division" as vector of strings.

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
			@param phase1 bool false indicates negation on in1
			@param phase2 bool false indicates negation on in2
			@param quotient std::vector<Monom> actually not used in current implementation
			@param quotientStrVec std::vector<std::string>
		*/
		void replaceANDDependingOnNegationsWithQuotients(varIndex replace, varIndex in1, varIndex in2, bool phase1, bool phase2, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec);
		
		/** Replace variable "replace" by the polynomial for a logic AND gate and additionally save the quotient of the "polynomial division" as vector of strings.   

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
			@param quotient std::vector<Monom> actually not used in current implementation
			@param quotientStrVec std::vector<std::string>
		*/	
		void replaceANDWithQuotients(varIndex replace, varIndex in1, varIndex in2, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec);
		
		/** Replace variable "replace" by the polynomial for a logic AND gate with inputs in1 and in2 where in1 is negated
			and additionally save the quotient of the "polynomial division" as vector of strings.         

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
			@param quotient std::vector<Monom> actually not used in current implementation
			@param quotientStrVec std::vector<std::string>
		*/
		void replaceANDOneNegationWithQuotients(varIndex replace, varIndex in1, varIndex in2, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec);
		
		/** Replace variable "replace" by the polynomial for a logic AND gate with both inputs negated 
			and additionally save the quotient of the "polynomial division" as vector of strings.      

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param in2 varIndex
			@param quotient std::vector<Monom> actually not used in current implementation
			@param quotientStrVec std::vector<std::string>
		*/
		void replaceANDDoubleNegationWithQuotients(varIndex replace, varIndex in1, varIndex in2, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec);
		
		/** Replace variable "replace" by the polynomial for a logic NOT gate with input in1
			and additionally save the quotient of the "polynomial division" as vector of strings.  

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param quotient std::vector<Monom> actually not used in current implementation
			@param quotientStrVec std::vector<std::string>
		*/
		void replaceNOTWithQuotients(varIndex replace, varIndex in1, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec);
		
		/** Replace variable "replace" by the polynomial for a logic BUFFER gate with input in1
			and additionally save the quotient of the "polynomial division" as vector of strings.  

			@param replace varIndex variable to replace
			@param in1 varIndex
			@param quotient std::vector<Monom> actually not used in current implementation
			@param quotientStrVec std::vector<std::string>
		*/
		void replaceBUFFERWithQuotients(varIndex replace, varIndex in1, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec);
		
		//************************ Getters and Setters.  ***************************************//
		
		/** Get a pointer to the set of monomials.  

			@return std::set<Monom>*
		*/
		const std::set<Monom>* getSet() const;
		
		/** Get a pointer to the beginning of RefList.  

			@return MyList*
		*/
		MyList* getRefList();
		
		/** Get a pointer to the phases vector.  

			@return std::vector<bool>*
		*/
		std::vector<bool>* getPhases();
		
		/** Set the phases vector.  

			@param std::vector<bool>
		*/
		void setPhases(std::vector<bool>& newPhases);
		
		/** Get varSize of the polynomial which indicates the max. variable index which can be saved in the polynomial. 

			@return size_t
		*/
		size_t getVarSize();
		
		/** Get the size (number of monomials) of the polynomial.

			@return size_t
		*/
		size_t size();
		
		//************************ Functions for printing out polynomial or monomials.  ***************************************//
		
		/** Print polynomial to standard output.

			@param stdout std::ostream
			@param obj Polynom
			@return std::ostream
		*/
		friend std::ostream& operator<<(std::ostream& stdout, const Polynom& obj);
		
		/** Return string of the polynomial.

			@return std::string
		*/
		const std::string to_string() const;
		
		/** Return string of the polynomial.

			@return std::string
		*/
		std::string to_string_opt() const;
		
		/** Return string of the polynomial, in reverse order.

			@return std::string
		*/
		const std::string to_string_reverse() const;
		
		/** Return string of the polynomial, indicating positive phases of variables by "x" and negative phases of variables by "f".

			@return std::string
		*/
		std::string to_string_with_phases() const;
		
		/** Return string of the polynomial, indicating positive phases of variables by "x" and negative phases of variables by "f".

			@return std::string
		*/
		std::string to_string_with_phases_opt() const;
		
		/** Replace all occurences of variable var in string str by its negation.
			"x_i" gets replaced by "f_i".

			@param str std::string
			@param var varIndex
		*/
		static void replace_var_by_negation(std::string& str, varIndex var);
		
		/** Return string of monomial, indicating positive phases of variables by "x" and negative phases of variables by "f".

			@return std::string
		*/
		std::string monToStringWithPhases(Monom mon) const;
		
		/** Return string of monomial.

			@return std::string
		*/
		std::string monToStringOpt(Monom mon) const;
		
		/** Return string of monomial, indicating positive phases of variables by "x" and negative phases of variables by "f".

			@return std::string
		*/
		std::string monToStringWithPhasesOpt(Monom mon) const;
		
		//*********************** Other helper functions.  ******************************************************//
		
		/** Resize the polynomials variable range. 
			Warning: Polynomial contents gets deleted and an empty polynomial with the set variable range is created.

			@param std::string
		*/
		void resize(size_t varSize);
		
		/** Multiply two polynomials: p1*p2. 

			@param p1 Polynom
			@param p2 Polynom
			@return Polynom
		*/
		static Polynom multiplyPoly(Polynom& p1, Polynom& p2);

		/** Find shortest monomial in the polynomial. 
			This monomial represents the shortest model for which the polynomial evaluates to a value unequal 0.
			If multiple shortest monomials exist, returns only one of them.

			@return Monom
		*/
		Monom getShortestModel();
		
		//*********************** Functions for controlling proof generation.  ******************************************************//
		
		/** Returns strings for the current polynomial and the modulo coefficient used for modulo reduction (if set).  

			@param std::pair<std::string, std::string> first string is current polynomial, second string is modulo reduction coefficient
		*/
		std::pair<std::string, std::string> writeOutStartingPoly();
		
		/** Set proof generation mode. If true, proof steps will be saved for substitution steps.  

			@param mode bool
		*/
		void setProofGenerationMode(bool mode);
		
		/** Start proof generation by setting the filenames where the proof will be written to.  

			@param polyFile std::string Filename where substitution polynomials are written to
			@param proofFile std::string Filename where the actual proof steps in PAC format is written to
		*/
		void startProofGeneration(std::string polyFile, std::string proofFile);	
		
		/** Helper function to convert substitution step into a PAC proof axiom.  

			@param replace varIndex
			@param mons std::list<Monom>
		*/
		std::string writeReplacementAxiom(varIndex replace, std::list<Monom>& mons);
		
		/** Helper function to convert substitution step into a PAC proof axiom.  

			@param replace varIndex
			@param mons std::set<Monom>
		*/
		std::string writeReplacementAxiom(varIndex replace, std::set<Monom>* mons);

	private:
		// Polynomial consists of two data structures: 1) Set of all monomials  2) List of all reference to monomials for every variable.
		std::set<Monom> polySet;
		MyList* refList;
		
		// Helping variable for remembering the variable range of polynomial. Only variables until varSize can be saved in the polynomial.
		size_t varSize;

		// Phase vector. Save which variables are currently in negated phase. False: Variable is negated. True: Variable is not negated.
		std::vector<bool> phases;

		// Mod reduction helpers.
		bool modReductionEnabled = false;
		mpz_class coefModReduction = 0;
		
		// Activating/deactivating proof writing.
		bool proofEnabled = false;
};

#endif /* POLYNOM_H_ */
