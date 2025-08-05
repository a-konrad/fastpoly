/*------------------------------------------------------------------------*/
/*! \file proof_writer.cpp
    \brief contains helper functions for generating PAC proofs.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

#include "proof_writer.h"
#include "polynom.h"

//****************************************************************************************/
// Global variables
int axiomNum = 0;
std::string polyfilename;
std::string prooffilename;
mpz_class modCoefProof = 0;
bool firstline = true;

//****************************************************************************************/
void set_proof_filenames(std::string polyname, std::string proofname) {
  	polyfilename = polyname;
  	prooffilename = proofname;
}

//****************************************************************************************/
void writeStartPolyToFile(std::pair<std::string, std::string> inputpair, int maxVarIndex) {
	axiomNum = 0;
	firstline = true;
	std::ofstream outputPolys;
	std::string outputPolysName = polyfilename;
	outputPolys.open(outputPolysName, std::ofstream::trunc); // Create a file for writing first axioms for the PAC proof.
	if (outputPolys.is_open()) {
		outputPolys << axiomNum++ << " " << maxVarIndex << ";" << std::endl;
		outputPolys << axiomNum++ << " " << inputpair.second << ";" << std::endl;
		outputPolys << axiomNum++ << " " << convertPolyStringToPACFormat(inputpair.first) << ";" << std::endl;
	}
	std::ofstream outputProof;
	std::string outputProofName = prooffilename;
	outputProof.open(outputProofName, std::ofstream::trunc);  // Open the file once so contents gets deleted if it already exists.
}

//****************************************************************************************/
void writeNewPolyAxiom(std::string axiomStr) {
	std::ofstream outputPolys;
	std::string outputPolysName = polyfilename;
	outputPolys.open(outputPolysName, std::ofstream::app); // Open file for appending axioms for the PAC proof.
	if (outputPolys.is_open()) {
		outputPolys << axiomNum++ << " " << convertPolyStringToPACFormat(axiomStr) << ";" << std::endl;
	}
}

//***************************************************************************************
std::string convertPolyStringToPACFormat(std::string subStr) {
	subStr.erase(std::remove(subStr.begin(), subStr.end(), '['), subStr.end());
	subStr.erase(std::remove(subStr.begin(), subStr.end(), ']'), subStr.end());
	subStr.erase(std::remove(subStr.begin(), subStr.end(), ' '), subStr.end());
	size_t pos = 0;
	while ((pos = subStr.find("+-")) != std::string::npos) {
        subStr.replace(pos, 2, "-");
    }
	while ((pos = subStr.find("*-")) != std::string::npos) {
        subStr.replace(pos, 2, "-");
    }
	while ((pos = subStr.find("*+")) != std::string::npos) {
        subStr.replace(pos, 2, "+");
    }
    while ((pos = subStr.find("-x")) != std::string::npos) {
        subStr.replace(pos, 2, "-1*x");
    }
	if (subStr == "1*") subStr.erase(1,1);
	return subStr;
}

//***************************************************************************************
void writePolysIntoPACProof(std::string inputName, std::string outputName) {
	Polynom pol;
	init_spec_from_PAC(pol, inputName);
	reduce_poly_with_proof(pol, inputName, outputName);
}

//****************************************************************************************/
void init_spec_from_PAC(Polynom & spec, std::string filename) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cout << "Error opening file " << filename << ". Make sure the filename is correct." << std::endl;
  }
  int lineNum = 0;
  unsigned maxVarNum = 0;
  std::string line;
  mpz_class modCoef = 0;
  while (std::getline(infile, line)) {
    ++lineNum;
    removeLineNumAndSemicolon(line);
    if (lineNum == 1) {
      maxVarNum = std::stoi(line);
      spec.resize(maxVarNum + 1);
    } else if (lineNum == 2) {
	  modCoefProof = mpz_class(line);
    } else if (lineNum == 3) {
      // Create spec poly.
      read_spec_poly_from_PAC(spec, line);
    } else {
     break;
    }
  }
  infile.close();
}

//****************************************************************************************/
void read_spec_poly_from_PAC(Polynom & spec, std::string line) {
  std::vector<std::string> monomialStrings;
  std::vector<std::string> variableStrings;
  std::regex monDelimit("([+-]?[^+-]+)");
  std::sregex_token_iterator it(line.begin(), line.end(), monDelimit);
  std::sregex_token_iterator end;
  for (; it != end; ++it) {
    monomialStrings.push_back(it->str());
  }
  for (auto& elem: monomialStrings) {
    variableStrings.clear();
    std::regex varDelimit("[*]");
    std::sregex_token_iterator it(elem.begin(), elem.end(), varDelimit, -1);
    std::sregex_token_iterator end;
    for (; it != end; ++it) {
      variableStrings.push_back(it->str());
    }
    std::string currStr;
    mpz_class tmpCoef;
	int* tmpVars = new int[variableStrings.size() - 1];
	int tmpSum = 0;
	for (size_t i=0; i < variableStrings.size(); ++i) {
		currStr = variableStrings[i];
		if (i==0) {
			size_t pos = 0;
			if ((pos = currStr.find("+")) != std::string::npos) currStr.erase(pos, 1);
			tmpCoef = mpz_class(currStr);
		} else {
			currStr.erase(currStr.find("x"), 1);
			tmpVars[i-1] = std::stoi(currStr);
			tmpSum += tmpVars[i-1];
		}
	}
	std::sort(tmpVars, tmpVars + variableStrings.size() - 1);
	Monom tmpMon(tmpVars, variableStrings.size() - 1, tmpSum, tmpCoef);
	delete[] tmpVars;
	spec.addMonom(tmpMon);
  }
}

//****************************************************************************************/
void reduce_poly_with_proof(Polynom& spec, std::string inputname, std::string outputname) {
  std::ifstream infile(inputname);
  if (!infile.is_open()) {
    std::cout << "Error opening file " << inputname << ". Make sure the filename is correct." << std::endl;
  }
  int lineNum = 0;
  unsigned maxVarNum = 0;
  std::string line;
  unsigned maxSize = 0;
  while (std::getline(infile, line)) {
    ++lineNum;
    if (lineNum < 4) continue;
    removeLineNumAndSemicolon(line);
    reduce_by_one_line_with_proof(spec, line, outputname, lineNum);
    if (maxSize < spec.size()) maxSize = spec.size();
    std::cout << "Current step: " << lineNum - 3 << " with poly.size: " << spec.size() << std::endl;
  }
  std::cout << "Steps completed." << std::endl;
  std::cout << "Max. Size was " << maxSize << std::endl;
}

//****************************************************************************************/
void reduce_by_one_line_with_proof(Polynom & spec, std::string line, std::string outputname, int lineNum) {
  std::vector<std::string> monomialStrings;
  std::vector<std::string> variableStrings;
  std::regex monDelimit("([+-]?[^+-]+)");
  std::sregex_token_iterator it(line.begin(), line.end(), monDelimit);
  std::sregex_token_iterator end;
  for (; it != end; ++it) {
    monomialStrings.push_back(it->str());
  }
  int leadingVar = 0;
  bool firstMonomial = true;
  std::list<Monom> tail;
  for (int j=0; j < monomialStrings.size(); ++j) {
    variableStrings.clear();
    std::regex varDelimit("[*]");
    std::sregex_token_iterator it(monomialStrings[j].begin(), monomialStrings[j].end(), varDelimit, -1);
    std::sregex_token_iterator end;
    for (; it != end; ++it) {
      variableStrings.push_back(it->str());
    }
    std::string currStr;
    mpz_class tmpCoef;
	int* tmpVars = new int[variableStrings.size() - 1];
	int tmpSum = 0;
	for (size_t i=0; i < variableStrings.size(); ++i) {
		currStr = variableStrings[i];
		if (i==0) {  // coefficient
			size_t pos = 0;
			if ((pos = currStr.find("+")) != std::string::npos) currStr.erase(pos, 1);
			tmpCoef = mpz_class(currStr);
		} else {  // variables
			currStr.erase(currStr.find("x"), 1);
			tmpVars[i-1] = std::stoi(currStr);
			tmpSum += tmpVars[i-1];
			if (firstMonomial) leadingVar = tmpVars[i-1];
		}
	}
	Monom tmpMon(tmpVars, variableStrings.size() - 1, tmpSum, tmpCoef);
	delete[] tmpVars;
	if (!firstMonomial) tail.push_back(tmpMon);
	firstMonomial = false;
  }
  std::vector<std::string> quotientStrVec;
  std::vector<Monom> quotient;
  spec.replaceVarWithQuotients(leadingVar, tail, quotient , quotientStrVec);
  writeOneLineIntoProof(outputname, lineNum - 1, spec, quotientStrVec);
}

/****************************************************************************************/
void writeOneLineIntoProof(std::string outputname, int usedAxiom, Polynom& pol, std::vector<std::string>& quotientStrVec) {
	std::string result;
	std::string quotientStr;
	for (size_t i=0; i < quotientStrVec.size(); ++i) {
		quotientStr.append(quotientStrVec[i]);
		if ((i != quotientStrVec.size() - 1) && (quotientStrVec[i+1][0] != '-')) quotientStr.append("+");
	}
	result.append(std::to_string(axiomNum++));
	result.append(" % ");
	result.append(std::to_string(usedAxiom));
	result.append(" *(");
	result.append(quotientStr);  //result.append(returnedQuotient);
	result.append(") + ");
	bool writeDelete = true;
	if (!firstline) result.append(std::to_string(axiomNum - 2));
	else { 
		result.append("2");
		firstline = false;
		writeDelete = false;
	}
	result.append(", ");
	result.append(pol.to_string_opt());  // result.append(returnedRemainder);
	result.append(";");
	if (modCoefProof > 0) result = addModReductionStep(pol, result, modCoefProof);
	// Create output stream and write axiom into it. Also write deletion of last axiom since it will not be used anymore.
	std::ofstream outputProof;
	std::string outputProofName = outputname;
	outputProof.open(outputProofName, std::ofstream::app); // Create a file for writing first axioms for the PAC proof.
	if (outputProof.is_open()) {
		outputProof << result << std::endl;
		if (writeDelete) outputProof << std::to_string(axiomNum - 2) << " d;" << std::endl;
	}
}

//***************************************************************************************
std::string addModReductionStep(Polynom& poly, std::string currStr, mpz_class modReduction) {
	std::string result = currStr;
	std::string multAxiom = "1";
	std::string returnedQuotient, addedAxiom, returnedRemainder;
	std::vector<std::string> quotientStrVec = poly.modReductionWithQuotientStr(modReduction);
	if (quotientStrVec.empty()) return result;
	for (size_t i=0; i < quotientStrVec.size(); ++i) {
		returnedQuotient.append(quotientStrVec[i]);
		if ((i != quotientStrVec.size() - 1) && (quotientStrVec[i+1][0] != '-')) returnedQuotient.append("+");
	}
	// Remove the current poly from currStr before adding mod reduction quotient and new resulting poly.
	result.erase(result.begin() + result.find(","), result.end());
	result.append(" + ");
	result.append(multAxiom);
	result.append(" *(");
	result.append(returnedQuotient);
	result.append(")");  // result.append(") + ");
	result.append(", ");
	result.append(poly.to_string_opt());  // result.append(returnedRemainder);
	result.append(";");
	return result;
}

//****************************************************************************************/
void removeLineNumAndSemicolon(std::string& line) {
    line.erase(std::remove(line.begin(), line.end(), ';'), line.end());
	line.erase(0, line.find(" ") + 1);
}

