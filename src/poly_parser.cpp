/*------------------------------------------------------------------------*/
/*! \file poly_parser.cpp
    \brief contains helper functions for reading polynomials and substitution
    steps from file input.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

#include "poly_parser.h"

//****************************************************************************************/
void init_spec(Polynom & spec, std::string filename) {
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
    if (lineNum == 1) {
      maxVarNum = std::stoi(line);
      spec.resize(maxVarNum + 1);
    } else if (lineNum == 2) {
      modCoef = mpz_class(line);
    } else if (lineNum == 3) {
      // Create spec poly.
      read_spec_poly(spec, line);
      if (modCoef > 0) {
        spec.setModReduction(true);
        spec.setModReductionNumber(modCoef);
        spec.modReducePoly(modCoef);
      }
    } else {
     break;
    }
  }
  infile.close();
}

//****************************************************************************************/
void read_spec_poly(Polynom & spec, std::string line) {
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
void reduce_poly(Polynom & spec, std::string filename) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cout << "Error opening file " << filename << ". Make sure the filename is correct." << std::endl;
  }
  int lineNum = 0;
  unsigned maxVarNum = 0;
  std::string line;
  unsigned maxSize = 0;
  while (std::getline(infile, line)) {
    ++lineNum;
    if (lineNum < 4) continue;
    reduce_by_one_line(spec, line);
    if (maxSize < spec.size()) maxSize = spec.size();
    std::cout << "Current step: " << lineNum - 3 << " with poly.size: " << spec.size() << std::endl;
  }
  std::cout << "Steps completed." << std::endl;
  std::cout << "Max. Size was " << maxSize << std::endl;
}

//****************************************************************************************/

void reduce_by_one_line(Polynom & spec, std::string line) {
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
  spec.replaceVar(leadingVar, tail);
}
