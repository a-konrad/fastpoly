/*------------------------------------------------------------------------*/
/*! \file monom.h
    \brief contains the class Monom for representing monomials.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

// std includes.
#include <string>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <sstream>

// Gnu multiprecision library.
#include <gmpxx.h>

// Local includes.
#include "mylist.h"

#ifndef MONOM_H_
#define MONOM_H_

// Integers used for variable indices. 
typedef int varIndex;

// Class to represent a monomial.
class Monom {

	friend class Polynom;
	public:	
		//*********************** Constructors  ******************************************************//
		
		/** Default Constructor. */
		Monom();
		
		/** Copy constructor.

			@param old Monom to copy
		*/
		Monom(const Monom& old);  
		
		/** Assignment operator.

			@param old Monom to copy
		*/
		Monom& operator = (const Monom& old);
		
		/** Constructor for monomial with one variable and coefficient 1. 

			@param index varIndex
		*/
		Monom(varIndex index);
		
		/** Constructor for monomial with two variables and coefficient 1. 

			@param index1 varIndex
			@param index2 varIndex
		*/
		Monom(varIndex index1, varIndex index2);
		
		/** Constructor for monomial from varIndex array. This constructor sorts the varIndex array and removes duplicates.

			@param myints[] varIndex
			@param size integer size of array argument
		*/
		Monom(varIndex myints[], int size);
		
		/** Constructor used internally. This constructor expects the array to be already sorted and duplicate free.

			@param myints[] varIndex
			@param size integer size of array argument
			@param sum integer sum of all array elements
			@param mpz_class factor coefficient of monomial
		*/
		Monom(varIndex myints[], int size, int sum, mpz_class factor);
		
		/** Destructor. */
		virtual ~Monom();
		
		//*********************** Getters and Setters  ******************************************************//
		
		/** Getter for variables of monomial.

			@return array of varIndex variables
		*/
		varIndex* getVars() const;
		
		/** Getter for MyList::ListElement array. 

			@return array of list elements
		*/
		MyList::ListElement** getPtrs() const;
		
		/** Setter for MyList::ListElement array. 

			@param ptrs MyList::ListElement**
		*/
		void setPtrs(MyList::ListElement** ptrs);
		
		/** Getter for monomial size. 

			@return integer
		*/
		int getSize() const;
		
		/** Setter for monomial size. 

			@param size int
		*/
		void setSize(int size);
		
		/** Getter for monomial coefficient.

			@return mpz_class
		*/
		mpz_class getFactor() const;
		
		/** Setter for monomial coefficient. 

			@param fact mpz_class
		*/
		void setFactor(mpz_class fact);
		
		/** Getter for monomials variable sum.

			@return int
		*/
		int getSum() const;
		
		/** Setter for monomials variable sum. 

			@param newSum int
		*/
		void setSum(int newSum);
		
		//*********************** Comparison operators.  ******************************************************//
		
		bool operator>(const Monom &m1) const;
		bool operator<(const Monom &m1) const;
		bool operator==(const Monom &m1) const;
		bool operator!=(const Monom &m1) const;
		
		//*********************** Merge function.  ******************************************************//
		
		// Replacing one variable in monomial by a whole new monomial mon. Return resulting monomial.
		/** Key function for backward rewriting operations. Replaces a variable with a given monomial. 

			@param replace varIndex to be replaced
			@return mon Monom newly created
		*/
		Monom merge(varIndex replace, Monom mon);
		
		//*********************** Other helper functions.  ******************************************************//
		
		/** Calculate sum by adding up all variable indices.

			@return int
		*/
		int calculateSum() const;
		
		/** Return whether variable v is contained in the monomial.

			@param v varIndex
			@return bool
		*/
		bool containsVar(varIndex v);
		
		/** Multiply two monomials and return the resulting monomial.

			@param mon1 Monom
			@param mon2 Monom
			@return Monom
		*/
		static Monom multiply(Monom mon1, Monom mon2);
		
		/** Print a monomial to standard output.

			@param stdout std::ostream
			@param obj Monom
			@return std::ostream
		*/
		friend std::ostream& operator<<(std::ostream& stdout, const Monom& obj);
	
		/** Return string of the monomial.

			@return std::string
		*/
		const std::string to_string() const;
		
		/** Return string of the monomial backwards.

			@return std::string
		*/
		const std::string to_string_reverse() const;
	
	private:
		// Array for all variables of monomial. 
		varIndex *vars;
		
		// Helper variables to track size (# of variables) and sum (sum of indices) of monomial.
		int size;
		int sum;
		
		// Coefficient of monomial. GMP library used because coefficients fastly exceed int64 range.
		mutable mpz_class factor;
		
		// Pointer back to ListElement entry, used for enabling constant deletion of elements from the list.
		MyList::ListElement** ptrs;



};

#endif /* MONOM_H_ */
