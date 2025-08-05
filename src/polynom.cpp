/*------------------------------------------------------------------------*/
/*! \file polynom.cpp
    \brief contains the class Polynom for representing polynomials and 
    defining all polynomial manipulation functions.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

#include "polynom.h"

//***************************************************************************************
Polynom::Polynom(){
	this->refList = new MyList[1];
	this->varSize = 1; 
	this->phases = {true};
}

//***************************************************************************************
Polynom::Polynom(int varSize){
	this->refList = new MyList[varSize+1];
	this->varSize = varSize;
	this->phases = std::vector<bool>(varSize+1, true);
}

//***************************************************************************************
Polynom::Polynom(const Polynom& old){  // Copy constructor
	this->varSize = old.varSize;
	this->refList = new MyList[old.varSize];
	this->polySet.clear(); 
	for (std::set<Monom>::iterator it = old.polySet.begin(); it != old.polySet.end(); ++it) {
		this->addMonom(*it);
	}
	this->phases = old.phases;
}

//***************************************************************************************
Polynom& Polynom::operator=(const Polynom& other) {  // Assignment operator.
	if (this != &other) {
		this->varSize = other.varSize;
		delete[] this->refList;
		this->refList = new MyList[other.varSize];
		this->polySet.clear();
		for (std::set<Monom>::iterator it = other.polySet.begin(); it != other.polySet.end(); ++it) {
			this->addMonom(*it);
		}
		this->phases = other.phases;
	}
	return *this;
}

//***************************************************************************************
Polynom::~Polynom(){
	delete[] this->refList;
}

//***************************************************************************************
bool Polynom::addPolynom(const Polynom& other) {
	if (this->varSize < other.varSize) {
		std::cout << "Cant add big polynom to small polynom (considering variable range). " << std::endl;
		return false;
	} else {
		for (std::set<Monom>::iterator it = other.polySet.begin(); it != other.polySet.end(); ++it) {
				this->addMonom(*it);
		}
		return true;
	}
}

//***************************************************************************************
Monom* Polynom::addMonom(Monom mon){
	std::pair<std::set<Monom>::iterator,bool> ret;
	ret = this->polySet.insert(mon);
	if (ret.second == false) {  //Monom alredy exists. Just add the factor. Check for 0 factor monoms.
		if ((ret.first->getFactor() + mon.getFactor()) == 0) {
			this->eraseMonom(*ret.first);  // Erase monom if factor is set to 0.
			return NULL;
		} else {
			ret.first->factor = ret.first->getFactor() + mon.getFactor();
		}
	} else { // New monom inserted.
		varIndex* vars = ret.first->getVars();
		int size = ret.first->getSize();
		for (int i = 0; i < size; i++) {  // Add reference to newly inserted monomials.
			this->addRefVar(const_cast<Monom&>(*ret.first), vars[i], i);
		}	
	}
	if (this->modReductionEnabled) {
		mpz_mod(ret.first->factor.get_mpz_t(), ret.first->factor.get_mpz_t(), this->coefModReduction.get_mpz_t());
		if (ret.first->getFactor() == 0) {
			this->eraseMonom(*ret.first); // Erase monom if factor after mod reduction is set to 0.
			return NULL;
		}
	}
	return &const_cast<Monom&>(*ret.first);
}

//***************************************************************************************
void Polynom::eraseMonom(Monom mon) {
	varIndex var = 0;
	for (int i = 0; i < mon.getSize(); i++) {
		var = (mon.getVars())[i];
		this->refList[var].deleteElement(mon.getPtrs()[i]);
		
	}
	size_t deletedElements = 0;	
	deletedElements = this->polySet.erase(mon);
}

//***************************************************************************************
void Polynom::createMonom(varIndex index1, mpz_class coef) {
	Monom tmp(index1);
	if (coef != 0) tmp.setFactor(coef);
	this->addMonom(tmp);
}

//***************************************************************************************
void Polynom::createMonom(varIndex index1, varIndex index2, mpz_class coef) {
	Monom tmp;
	if (index2==0) tmp = Monom(index1);
	else tmp = Monom(index1, index2);
	if (coef != 0) tmp.setFactor(coef);
	this->addMonom(tmp);
}

//***************************************************************************************
void Polynom::createMonom(varIndex myints[], int size, mpz_class coef) {
	Monom tmp(myints, size);
	if (coef != 0) tmp.setFactor(coef);
	this->addMonom(tmp);
}

//***************************************************************************************
void Polynom::replaceANDDependingOnNegations(varIndex replace, varIndex in1, varIndex in2, bool phase1, bool phase2) {
	if (phase1) {
		if (phase2) {
			replaceAND(replace, in1, in2);
		} else {
			replaceANDOneNegation(replace, in2, in1);
		}
	} else if (phase2) {
			replaceANDOneNegation(replace, in1, in2);
		} else {
			replaceANDDoubleNegation(replace, in2, in1);
		}
}

//***************************************************************************************
void Polynom::replaceAND(varIndex replace, varIndex in1, varIndex in2) {
	varIndex tmp = -1;
	if (in1 == in2) {
		Monom oneMonom(in1);
		std::list<Monom> mons;
		mons.push_back(oneMonom);
		this->replaceVar(replace, mons);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
		}	
		Monom andMon(in1, in2);
		std::list<Monom> mons;
		mons.push_back(andMon);
		this->replaceVar(replace, mons);
	}
}

//***************************************************************************************
void Polynom::replaceANDOneNegation(varIndex replace, varIndex in1, varIndex in2) {
	// It is important, that in1 is the negated signal.
	varIndex tmp = -1;
	bool switched = false;
	if (in1 == in2) {
		Monom zeroMonom;
		zeroMonom.setFactor(0);
		std::list<Monom> mons;
		mons.push_back(zeroMonom);
		this->replaceVar(replace, mons);
	} else {
		if (in1 > in2) { // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
			switched = true;
		}
		Monom andMon(in1, in2);
		andMon.setFactor(-1);
		Monom andMon2;
		if (switched) andMon2 = Monom(in1);  // If signals were swapped, in1 is now the non inverted signal.
		else andMon2 = Monom(in2);
		std::list<Monom> mons;
		mons.push_back(andMon);
		mons.push_back(andMon2);
		this->replaceVar(replace, mons);
	}
}

//***************************************************************************************
void Polynom::replaceANDDoubleNegation(varIndex replace, varIndex in1, varIndex in2) {
	varIndex tmp = -1;
	if (in1 == in2) {
		replaceNOT(replace, in1);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
		}
		Monom andMon(in1, in2);
		Monom andMon2(in1);
		andMon2.setFactor(-1);
		Monom andMon3(in2);
		andMon3.setFactor(-1);
		Monom andMon4;
		andMon4.setFactor(1);
		std::list<Monom> mons;
		mons.push_back(andMon);
		mons.push_back(andMon2);
		mons.push_back(andMon3);
		mons.push_back(andMon4);
		this->replaceVar(replace, mons);
	}
}

//***************************************************************************************
void Polynom::replaceOR(varIndex replace, varIndex in1, varIndex in2) {
	varIndex tmp = -1;
	if (in1 == in2) {
		Monom oneMonom(in1);
		std::list<Monom> mons;
		mons.push_back(oneMonom);
		this->replaceVar(replace, mons);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;		
		}	
		Monom orMon1(in1);
		Monom orMon2(in2);
		Monom orMon3(in1, in2);
		orMon3.setFactor(-1);
		std::list<Monom> mons;
		mons.push_back(orMon1);
		mons.push_back(orMon2);
		mons.push_back(orMon3);
		this->replaceVar(replace, mons);
	}
}

//***************************************************************************************
void Polynom::replaceOROneNegation(varIndex replace, varIndex in1, varIndex in2) {
	// It is important, that in1 is the negation signal.
	varIndex tmp = -1;
	bool switched = false;
	if (in1 == in2) {
		Monom oneMonom;
		oneMonom.setFactor(1);
		std::list<Monom> mons;
		mons.push_back(oneMonom);
		this->replaceVar(replace, mons);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
			switched = true;
		}
		Monom orMon1;
		orMon1.setFactor(1);
		Monom orMon2;
		if (switched) orMon2 = Monom(in2);  // If switched, in2 is now the negation signal.
		else orMon2 = Monom(in1);
		orMon2.setFactor(-1);
		Monom orMon3(in1, in2);
		std::list<Monom> mons;
		mons.push_back(orMon1);
		mons.push_back(orMon2);
		mons.push_back(orMon3);
		this->replaceVar(replace, mons);
	}
}

//***************************************************************************************
void Polynom::replaceORDoubleNegation(varIndex replace, varIndex in1, varIndex in2) { 
	varIndex tmp = -1;
	if (in1 == in2) {
		replaceNOT(replace, in1);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
		}
		Monom orMon1;
		orMon1.setFactor(1);
		Monom orMon2(in1, in2);
		orMon2.setFactor(-1);
		std::list<Monom> mons;
		mons.push_back(orMon1);
		mons.push_back(orMon2);
		this->replaceVar(replace, mons);
	}
}

//***************************************************************************************
void Polynom::replaceXOR(varIndex replace, varIndex in1, varIndex in2) {
	varIndex tmp = -1;
	if (in1 == in2) {
		Monom empty;
		empty.setFactor(0);
		std::list<Monom> mons;
		mons.push_back(empty);
		this->replaceVar(replace, mons);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
		}	
		Monom xorMon1(in1);
		Monom xorMon2(in2);
		Monom xorMon3(in1, in2);
		xorMon3.setFactor(-2);
		std::list<Monom> mons;
		mons.push_back(xorMon1);
		mons.push_back(xorMon2);
		mons.push_back(xorMon3);
		this->replaceVar(replace, mons);
	}
}

//***************************************************************************************
void Polynom::replaceXOROneNegation(varIndex replace, varIndex in1, varIndex in2) {
	varIndex tmp = -1;
	if (in1 == in2) {
		Monom oneMonom;
		oneMonom.setFactor(1);
		std::list<Monom> mons;
		mons.push_back(oneMonom);
		this->replaceVar(replace, mons);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
		}
		Monom xorMon1(in1);
		xorMon1.setFactor(-1);
		Monom xorMon2(in2);
		xorMon2.setFactor(-1);
		Monom xorMon3(in1, in2);
		xorMon3.setFactor(2);
		Monom xorMon4;
		xorMon4.setFactor(1);
		std::list<Monom> mons;
		mons.push_back(xorMon1);
		mons.push_back(xorMon2);
		mons.push_back(xorMon3);
		mons.push_back(xorMon4);
		this->replaceVar(replace, mons);
	}
}

//***************************************************************************************
void Polynom::replaceNOT(varIndex replace, varIndex in1) {
	Monom notMon1(in1);
	notMon1.setFactor(-1);
	Monom notMon2;
	notMon2.setFactor(1);
	std::list<Monom> mons;
	mons.push_back(notMon1);
	mons.push_back(notMon2);
	this->replaceVar(replace, mons);
}

//***************************************************************************************
void Polynom::replaceBUFFER(varIndex replace, varIndex in1) {
	if (replace == in1) {std::cout << "replaceBuffer with same variables" << std::endl; return; }
	Monom notMon1(in1);
	notMon1.setFactor(1);
	std::list<Monom> mons;
	mons.push_back(notMon1);
	this->replaceVar(replace, mons);
}

//***************************************************************************************
void Polynom::replaceCON0(varIndex replace) {
	Monom notMon1;
	notMon1.setFactor(0);
	std::list<Monom> mons;
	mons.push_back(notMon1);
	this->replaceVar(replace, mons);
}

//***************************************************************************************
void Polynom::replaceCON1(varIndex replace) {
	Monom notMon1;
	notMon1.setFactor(1);
	std::list<Monom> mons;
	mons.push_back(notMon1);
	this->replaceVar(replace, mons);
}

//***************************************************************************************
void Polynom::replaceVarByPoly(varIndex replace, Polynom& poly) {
	this->replaceVar(replace, const_cast<std::set<Monom>*>(poly.getSet()));
}

//***************************************************************************************
void Polynom::replaceVar(varIndex replace, std::list<Monom>& mons) {
	if (this->proofEnabled) writeNewPolyAxiom(writeReplacementAxiom(replace, mons));
	Monom newMon;
	Monom oldMon;
	MyList::ListElement* nextElement;
	Monom* newMonPointer = NULL;
	std::pair<std::set<Monom>::iterator, bool> retPair;
	for (MyList::Iterator it=this->refList[replace].begin(); it != this->refList[replace].end(); it = nextElement) {
		if (refList[replace].isEmpty()) {
			std::cout << "Reflist Empty. Something went wrong." << std::endl;
			return;
		}
		oldMon = *(it.returnData());
		nextElement = it.returnElement()->next;  // Save address of next element before deleting current element.
		this->eraseMonom(*(it.returnData()));  // ATTENTION: Iterator on the current element gets invalid.
		it = this->refList[replace].begin();
		for (std::list<Monom>::iterator it2=mons.begin(); it2 != mons.end(); ++it2) {
			newMon = oldMon.merge(replace, *it2);
			if (newMon.getFactor() == 0) continue; // Dont add monom with factor 0. Only caused by XOR with same inputs.
			newMonPointer = this->addMonom(newMon);
		}
	}
}

//***************************************************************************************
void Polynom::replaceVar(varIndex replace, std::set<Monom>* mons) {
	if (this->proofEnabled) writeNewPolyAxiom(writeReplacementAxiom(replace, mons));
	Monom newMon;
	Monom oldMon;
	MyList::ListElement* nextElement;
	Monom* newMonPointer = NULL;
	std::pair<std::set<Monom>::iterator, bool> retPair;
	for (MyList::Iterator it=this->refList[replace].begin(); it != this->refList[replace].end(); it = nextElement) {
		if (refList[replace].isEmpty()) {
			std::cout << "Reflist Empty. Something went wrong." << std::endl;
			return;
		}
		oldMon = *(it.returnData());
		nextElement = it.returnElement()->next;  // Save address of next element before deleting current element.
		this->eraseMonom(*(it.returnData()));  // ATTENTION: Iterator on the current element gets invalid.
		it = this->refList[replace].begin();
		for (std::set<Monom>::iterator it2=mons->begin(); it2 != mons->end(); ++it2) {
			newMon = oldMon.merge(replace, *it2);
			if (newMon.getFactor() == 0) continue; // Dont add monom with factor 0. Only caused by XOR with same inputs.
			newMonPointer = this->addMonom(newMon);
		}
	}
}

//***************************************************************************************
void Polynom::replaceVarWithQuotients(varIndex replace, std::list<Monom>& mons, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec) {
	Monom newMon;
	Monom oldMon;
	Monom con1;
	con1.setFactor(1);
	MyList::ListElement* nextElement;
	Monom* newMonPointer = NULL;
	std::pair<std::set<Monom>::iterator, bool> retPair;
	for (MyList::Iterator it=this->refList[replace].begin(); it != this->refList[replace].end(); it = nextElement) {
		if (refList[replace].isEmpty()) {
			std::cout << "Reflist Empty. Something went wrong." << std::endl;
			return;
		}
		oldMon = *(it.returnData());
		nextElement = it.returnElement()->next;  // Save address of next element before deleting current element.
		this->eraseMonom(*(it.returnData()));  // ATTENTION: Iterator on the current element gets invalid.
		it = this->refList[replace].begin();
		quotientStrVec.push_back(this->monToStringOpt(oldMon.merge(replace, con1)));
		for (std::list<Monom>::iterator it2=mons.begin(); it2 != mons.end(); ++it2) {
			newMon = oldMon.merge(replace, *it2);
			if (newMon.getFactor() == 0) continue; // Dont add monom with factor 0. Only caused by XOR with same inputs.
			newMonPointer = this->addMonom(newMon);
		}
	}
}

//***************************************************************************************
void Polynom::replaceANDDependingOnNegationsWithQuotients(varIndex replace, varIndex in1, varIndex in2, bool phase1, bool phase2, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec) {
	if (phase1) {
		if (phase2) {
			replaceANDWithQuotients(replace, in1, in2, quotient, quotientStrVec);
		} else {
			replaceANDOneNegationWithQuotients(replace, in2, in1, quotient, quotientStrVec);
		}
	} else if (phase2) {
			replaceANDOneNegationWithQuotients(replace, in1, in2, quotient, quotientStrVec);
		} else {
			replaceANDDoubleNegationWithQuotients(replace, in2, in1, quotient, quotientStrVec);
		}
}

//***************************************************************************************
void Polynom::replaceANDWithQuotients(varIndex replace, varIndex in1, varIndex in2, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec) {
	varIndex tmp = -1;
	if (in1 == in2) {
		Monom oneMonom(in1);
		std::list<Monom> mons;
		mons.push_back(oneMonom);
		this->replaceVarWithQuotients(replace, mons, quotient, quotientStrVec);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
		}	
		Monom andMon(in1, in2);
		std::list<Monom> mons;
		mons.push_back(andMon);
		this->replaceVarWithQuotients(replace, mons, quotient, quotientStrVec);
	}
}

//***************************************************************************************
void Polynom::replaceANDOneNegationWithQuotients(varIndex replace, varIndex in1, varIndex in2, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec) {
	// It is important, that in1 is the negation signal.
	varIndex tmp = -1;
	bool switched = false;
	if (in1 == in2) {
		Monom zeroMonom;
		zeroMonom.setFactor(0);
		std::list<Monom> mons;
		mons.push_back(zeroMonom);
		this->replaceVarWithQuotients(replace, mons, quotient, quotientStrVec);
	} else {
		if (in1 > in2) { // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
			switched = true;
		}
		Monom andMon(in1, in2);
		andMon.setFactor(-1);
		Monom andMon2;
		if (switched) andMon2 = Monom(in1);  // If signals were swapped, in1 is now the not inverted signal.
		else andMon2 = Monom(in2);
		std::list<Monom> mons;
		mons.push_back(andMon);
		mons.push_back(andMon2);
		this->replaceVarWithQuotients(replace, mons, quotient, quotientStrVec);
	}
}

//***************************************************************************************
void Polynom::replaceANDDoubleNegationWithQuotients(varIndex replace, varIndex in1, varIndex in2, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec) {
	varIndex tmp = -1;
	if (in1 == in2) {
		replaceNOTWithQuotients(replace, in1, quotient, quotientStrVec);
	} else {
		if (in1 > in2) {  // Swap signals to assure in1 is the smaller one.
			tmp = in1;
			in1 = in2;
			in2 = tmp;
		}
		Monom andMon(in1, in2);
		Monom andMon2(in1);
		andMon2.setFactor(-1);
		Monom andMon3(in2);
		andMon3.setFactor(-1);
		Monom andMon4;
		andMon4.setFactor(1);
		std::list<Monom> mons;
		mons.push_back(andMon);
		mons.push_back(andMon2);
		mons.push_back(andMon3);
		mons.push_back(andMon4);
		this->replaceVarWithQuotients(replace, mons, quotient, quotientStrVec);
	}
}

//***************************************************************************************
void Polynom::replaceNOTWithQuotients(varIndex replace, varIndex in1, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec) {
	Monom notMon1(in1);
	notMon1.setFactor(-1);
	Monom notMon2;
	notMon2.setFactor(1);
	std::list<Monom> mons;
	mons.push_back(notMon1);
	mons.push_back(notMon2);
	this->replaceVarWithQuotients(replace, mons, quotient, quotientStrVec);
}

//***************************************************************************************
void Polynom::replaceBUFFERWithQuotients(varIndex replace, varIndex in1, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec) {
	if (replace == in1) {std::cout << "replaceBufferWithQuotients with same variables" << std::endl; return; }
	Monom notMon1(in1);
	notMon1.setFactor(1);
	std::list<Monom> mons;
	mons.push_back(notMon1);
	this->replaceVarWithQuotients(replace, mons, quotient, quotientStrVec);
}

//***************************************************************************************
void Polynom::negateVar(varIndex replace) {
	std::vector<Monom*> oldPointers = this->findContainingVar(replace);
	Monom tmpMon;
	Monom mergeMon;  // Create "1" monomial.
	mergeMon.setFactor(1);
	for (auto& elem: oldPointers) {
		tmpMon = *elem;  // Save old mon.
		this->eraseMonom(*elem);  // Erase old mon.
		tmpMon.setFactor(tmpMon.getFactor() * -1);  // Negate factor.
		this->addMonom(tmpMon);	// Insert mon with negated factor.
		tmpMon.setFactor(tmpMon.getFactor() * -1);  // Get previous factor back.
		this->addMonom(tmpMon.merge(replace, mergeMon));
	}
	this->phases[replace] = !this->phases[replace];
}

//***************************************************************************************
bool Polynom::testPhaseChangeSingleVariable(varIndex var) {
	size_t sizeBefore = this->size();
	bool phaseBefore = this->phases.at(var);
	if (containsVar(var)) {
		negateVar(var);
		if (sizeBefore <= this->size()) {
			this->phases.at(var) = phaseBefore;
			return false;
		} else {
			return true;
		}
	}
	return false;
}

//***************************************************************************************
bool Polynom::testPhaseChangeSingleVariableImproved(varIndex var) {
	size_t sizeBefore = this->size();
	std::vector<Monom*> oldPointers = this->findContainingVar(var);
	Monom tmpMon;  // Create constant "1" monomial.
	tmpMon.setFactor(1);
	std::vector<Monom> addedMons;
	addedMons.reserve(oldPointers.size());
	for (auto& elem: oldPointers) {
		addedMons.push_back(elem->merge(var, tmpMon));
		this->addMonom(addedMons.back());
		mpz_neg(elem->factor.get_mpz_t(), elem->factor.get_mpz_t());
	}
	if (sizeBefore <= this->size()) { // Revert the negation.
		for (auto& elem: addedMons) {
			mpz_neg(elem.factor.get_mpz_t(), elem.factor.get_mpz_t());
			this->addMonom(elem);  // Erase previosuly added monoms.
		}
		for (auto& elem: oldPointers) {
			mpz_neg(elem->factor.get_mpz_t(), elem->factor.get_mpz_t());
		}
		return false;
	} else {  // Keep the negation.
		this->phases[var] = !this->phases[var];
		return true;
	}
}

//***************************************************************************************
void Polynom::negateVarImproved(varIndex var) {
	size_t sizeBefore = this->size();
	std::vector<Monom*> oldPointers = this->findContainingVar(var);
	Monom tmpMon;  // Create "1" monomial.
	tmpMon.setFactor(1);
	for (auto& elem: oldPointers) {
		this->addMonom(elem->merge(var, tmpMon));
		mpz_neg(elem->factor.get_mpz_t(), elem->factor.get_mpz_t());
		if (modReductionEnabled) {  // If modulo reduction is enabled, apply modulo to the negated factors.
			mpz_mod(elem->factor.get_mpz_t(), elem->factor.get_mpz_t(), this->coefModReduction.get_mpz_t());
		}
	}
	this->phases[var] = !this->phases[var];
}

//***************************************************************************************
void Polynom::negateVarImprovedWithQuotient(varIndex var, std::vector<Monom>& quotient, std::vector<std::string>& quotientStrVec) {
	size_t sizeBefore = this->size();
	std::vector<Monom*> oldPointers = this->findContainingVar(var);
	Monom con1Mon;  // Create "1" monomial.
	con1Mon.setFactor(1);
	Monom tmpMon;
	for (auto& elem: oldPointers) {
		quotientStrVec.push_back(monToStringWithPhasesOpt(elem->merge(var, con1Mon)));
		this->addMonom(elem->merge(var, con1Mon));
		mpz_neg(elem->factor.get_mpz_t(), elem->factor.get_mpz_t());
	}
	this->phases[var] = !this->phases[var];
}

//***************************************************************************************
int Polynom::greedyPhaseChange() {
	int improvement = 0;
	size_t sizeStart = this->size();
	for (size_t i=0; i < getVarSize(); ++i) {
		testPhaseChangeSingleVariableImproved(i);
	}
	improvement = (sizeStart - this->size());
	return improvement;
}

//***************************************************************************************
int Polynom::greedyPhaseChangeBackward() {
	int improvement = 0;
	size_t sizeStart = this->size();
	for (int i = getVarSize() - 1; i >= 0; --i) {
		testPhaseChangeSingleVariableImproved(i);
	}
	improvement = (sizeStart - this->size());
	return improvement;
}

//***************************************************************************************
int Polynom::greedyPhaseChangeCustom(std::vector<varIndex>& signalsToChange) {
	int improvement = 0;
	size_t sizeStart = this->size();
	for (size_t i=0; i < signalsToChange.size(); ++i) {
		testPhaseChangeSingleVariableImproved(signalsToChange.at(i));
	}
	improvement = (sizeStart - this->size());
	return improvement;
}

//***************************************************************************************
int Polynom::greedyPhaseChangeCustom(std::list<varIndex>& signalsToChange) {
	int improvement = 0;
	size_t sizeStart = this->size();
	for (std::list<varIndex>::iterator it= signalsToChange.begin(); it != signalsToChange.end(); ++it) {
		testPhaseChangeSingleVariableImproved(*it);
	}
	improvement = (sizeStart - this->size());
	return improvement;
}

//***************************************************************************************
int Polynom::greedyPhaseChangeCustom(std::list<varIndex>& signalsToChange, std::list<uint32_t>& changedPhases) {
	int improvement = 0;
	size_t sizeStart = this->size();
	bool changed;
	for (std::list<varIndex>::iterator it= signalsToChange.begin(); it != signalsToChange.end(); ++it) {
		changed = testPhaseChangeSingleVariableImproved(*it);
		if (changed) {
			changedPhases.push_back(*it);
		}
	}
	improvement = (sizeStart - this->size());
	return improvement;
}

//***************************************************************************************
void Polynom::reportVarPhases() {
	std::cout << "Varialbes with face 0 are following: " << std::endl;
	for (size_t i=0; i < getVarSize(); ++i) {
		if (this->phases[i] == false) std::cout << "x" << i << " phase is " << this->phases[i] << std::endl;
	}
}

//****************************************************************************************************************************
std::vector<Monom*> Polynom::findContainingVar(varIndex var) {
	std::vector<Monom*> resultVec;
	int minListLength = INT_MAX;
	varIndex minListVar = -1;
	// Retrieve all monomials from refList which contain var.
	bool add;
	for (MyList::Iterator it=this->refList[var].begin(); it != this->refList[var].end(); it++) {
		resultVec.push_back(it.returnData());
	}
	return resultVec;
}

//****************************************************************************************************************************
std::vector<Monom*> Polynom::findContaining(Monom& mon) {
	std::vector<Monom*> resultVec;
	int minListLength = INT_MAX;
	varIndex minListVar = -1;
	for (size_t i=0; i < mon.getSize(); i++) {  // Get the variable with shortest refList. 
		if (mon.getVars()[i] > this->getVarSize()) {
			std::cout << "Error in findContaining(): Monomial to find includes a variable out ouf range of this polynomial variable range." << std::endl;
			return resultVec;
		}
		if ((getRefList()[mon.getVars()[i]]).getSize() < minListLength) {
			minListVar = mon.getVars()[i]; 
			minListLength = (getRefList()[mon.getVars()[i]]).getSize();
			if (minListLength == 0) return resultVec;  // If one refList length is zero, this variable is not contained, so mon cannot be contained in polynomial.
		}
	}
	// Retrieve all monomials from shortest refList which contain mon.
	bool add;
	for (MyList::Iterator it=this->refList[minListVar].begin(); it != this->refList[minListVar].end(); it++) {
		add = true;
		for (size_t i=0; i < mon.getSize(); i++) {
			if (mon.getVars()[i] == minListVar) continue;  // Not needed to check for minListVar, since it is the list of this variable.
			if (it.returnData()->containsVar(mon.getVars()[i]) == false) {
				add = false;
				break;
			}
		}
		if (add) resultVec.push_back(it.returnData());
	}
	return resultVec;
}


//****************************************************************************************************************************
Monom* Polynom::findExact(Monom& mon) {
	// First check if special case: mon is the empty monomial (only a coefficient without variables).
	Monom* temp;
	if (mon.getSize() == 0) {  // If size=0 it is the empty monomial which is always first in polySet
		temp = &const_cast<Monom&>(*this->polySet.begin());
		if (temp->getSize() == 0) return &const_cast<Monom&>(*this->polySet.begin());
		else return NULL;
	}
	int minListLength = INT_MAX;
	varIndex minListVar = -1;
	for (size_t i=0; i < mon.getSize(); i++) {  // Get the variable with shortest refList.  
		if (mon.getVars()[i] > this->getVarSize()) {
			std::cout << "Error in findExact(): Monomial to find includes a variable out ouf range of this polynomial variable range." << std::endl;
			return NULL;
		}
		if ((getRefList()[mon.getVars()[i]]).getSize() < minListLength) {
			minListVar = mon.getVars()[i]; 
			minListLength = (getRefList()[mon.getVars()[i]]).getSize();
			if (minListLength == 0) return NULL;  // If one refList length is zero, this variable is not contained, so mon cannot be contained in polynomial.
		}
	}
	// Find exact monomial mon from shortest refList.
	for (MyList::Iterator it=this->refList[minListVar].begin(); it != this->refList[minListVar].end(); it++) {
		if (*(it.returnData()) == mon) {
			return it.returnData();
		} 
	}
	return NULL;  // This case should never happen.
}

//****************************************************************************************************************************
int Polynom::phaseChangeEffectOnMonom(Monom& mon, varIndex var) {
	// First check if special case: mon is the empty monomial (only a coefficient without variables).
	Monom* temp;
	if (mon.getSize() == 0) {  // If size=0 it is the empty monomial which is always first in polySet
		return 0;
	}
	int minListLength = INT_MAX;
	varIndex minListVar = -1;
	for (size_t i=0; i < mon.getSize(); i++) {  // Get the variable with shortest refList.
		if (mon.getVars()[i] > this->getVarSize()) {
			std::cout << "Error in findExact(): Monomial to find includes a variable out ouf range of this polynomial variable range." << std::endl;
			return NULL;
		}
		if ((getRefList()[mon.getVars()[i]]).getSize() < minListLength) {
			minListVar = mon.getVars()[i];
			minListLength = (getRefList()[mon.getVars()[i]]).getSize();
			if (minListLength == 0) return NULL;  // If one refList length is zero, this variable is not contained, so mon cannot be contained in polynomial.
		}
	}
	// Find monomial which is the same except that var is missing. Use sum to faster find candidates.
	int64_t findSum = mon.getSum() - var;
	size_t findSize = mon.getSize() - 1;
	int polySizeChange = 1;  // If monomial not found we add 1 to the poly size. If found we either dont change size or reduce by 1.
	for (MyList::Iterator it=this->refList[minListVar].begin(); it != this->refList[minListVar].end(); it++) {
		if ((it.returnData())->getSum() != findSum) continue;
		if ((it.returnData())->getSize() != findSize) continue;
		varIndex currVar;
		varIndex monVar;
		int offset = 0;
		if (mon.getVars()[0] == var) ++offset;
		bool conOuter = false;
		for (size_t monPos=0; monPos < (it.returnData())->getSize(); ++monPos) {
			currVar = (it.returnData())->getVars()[monPos];
			if (offset == 0 & currVar > var) ++offset;
			monVar = mon.getVars()[monPos + offset];
			if (currVar != monVar) { conOuter = true; break; }  // This is not the searched monomial.
		}
		if (conOuter) continue;  // Check next monomial.
		// Searched monomial found. Poly size is not increased. If coefficient is negative of mon we even reduce poly size by 1.
		polySizeChange = 0;
		if ((it.returnData())->factor == (-1 * mon.factor)) polySizeChange = -1;
		if ((it.returnData())->factor + mon.factor == this->coefModReduction) polySizeChange = -1;
		break;
	}
	return polySizeChange;  // This case should never happen.
}


//****************************************************************************************************************************
bool Polynom::containsVar(varIndex var) {
	if (var > this->varSize) {
		std::cout << "Searched varIndex exceeds variable range of polynomial." << std::endl;
		return false;
	} else {
		if (this->refList[var].getSize() == 0) return false;
		else return true;
	}
}

//***************************************************************************************
void Polynom::addRefVar(Monom& mon, varIndex index, int i) {
	(mon.ptrs)[i] = this->refList[index].add(&mon);
	return;
}

//***************************************************************************************
const std::set<Monom>* Polynom::getSet() const {
	return &this->polySet;
}

//***************************************************************************************
MyList* Polynom::getRefList() {
	return this->refList;
}


//***************************************************************************************
std::vector<bool>* Polynom::getPhases() {
	return &this->phases;
}

//***************************************************************************************
void Polynom::setPhases(std::vector<bool>& newPhases) {
	this->phases = newPhases;
}

//***************************************************************************************
size_t Polynom::getVarSize() {
	return this->varSize;
}

//***************************************************************************************
size_t Polynom::size() {
	return this->polySet.size();
}


//***************************************************************************************
void Polynom::resize(size_t varSize) {
	delete[] this->refList;
	this->polySet.clear();
	this->refList = new MyList[varSize];
	this->phases= std::vector<bool>(varSize+1, true);
	this->varSize = varSize;
}

//***************************************************************************************
Polynom Polynom::multiplyPoly(Polynom& p1, Polynom& p2) {
	int maxSize = 0;
	if (p1.getVarSize() < p2.getVarSize()) {
		maxSize = p2.getVarSize();
	} else {
		maxSize = p1.getVarSize();
	}
	Polynom mult(maxSize);
	
	Monom temp;
	for (std::set<Monom>::iterator it= p1.polySet.begin(); it != p1.polySet.end(); ++it) {
		for (std::set<Monom>::iterator it2= p2.polySet.begin(); it2 != p2.polySet.end(); ++it2) {
			temp = Monom::multiply(*it, *it2);
			mult.addMonom(temp);
		}
	}
	
	return mult;
}

//***************************************************************************************
void Polynom::parsePolyFromString(std::string inputStr) {
	std::vector<std::string> monomialStrings;
	std::vector<std::string> variableStrings;
	std::regex monDelimit("([+-]?[^+-]+)");
	std::sregex_token_iterator it(inputStr.begin(), inputStr.end(), monDelimit);
	std::sregex_token_iterator end;
	for (; it != end; ++it) {
        monomialStrings.push_back(it->str());
    }
	int highestIndex = 0;
	std::vector<Monom> monomials;

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
		Monom tmpMon(tmpVars, variableStrings.size() - 1, tmpSum, tmpCoef);
		delete[] tmpVars;
		this->addMonom(tmpMon);
	}
}

//***************************************************************************************
std::ostream& operator<<(std::ostream& stdout, const Polynom& obj) {
    std::string start = "", end = "", delim = " + ";
    std::string s;
    int size = obj.getSet()->size();
    if (!obj.getSet()->empty()){
        s += start;
        int num = 0;
        for (std::set<Monom>::iterator it=obj.getSet()->begin(); it != obj.getSet()->end(); ++it) {
            s += (*it).to_string();
            num += 1;
            if (num != size) s += delim;
        }
        s += end;
        //Add factor at the end
    }
    else{
        s += "0 (empty polynomial)";
    }
    stdout << s;
    return stdout;
}

//***************************************************************************************
const std::string Polynom::to_string() const {
    std::ostringstream ss;
    ss << *this;
    return ss.str();
}

//***************************************************************************************
const std::string Polynom::to_string_reverse() const {
    std::string start = "", end = "", delim = " + ";
    std::string s;
    int size = this->getSet()->size();
    if (!this->getSet()->empty()){
        s += start;
        int num = 0;
        for (std::set<Monom>::reverse_iterator it=this->getSet()->rbegin(); it != this->getSet()->rend(); ++it) {
            s += (*it).to_string_reverse();
            num += 1;
            if (num != size) s += delim;
        }
        s += end;
        //Add factor at the end
    }
    else{
        s += "0 (empty polynomial)";
    }
    return s;
}

//***************************************************************************************
std::string Polynom::to_string_with_phases() const {
    std::string start = "", end = "", delim = " + ";
    std::string str;
	std::string tmp;
    int size = this->getSet()->size();
    if (!this->getSet()->empty()){
        str += start;
        int num = 0;
        for (std::set<Monom>::iterator it=this->getSet()->begin(); it != this->getSet()->end(); ++it) {

        	tmp = monToStringWithPhases(*it);
			str += tmp;
            num += 1;
            if (num != size) str += delim;
        }
        str += end;
        //Add factor at the end
    }
    else{
        str += "0";
    }
    return str;
}

//***************************************************************************************
std::string Polynom::monToStringWithPhases(Monom mon) const {
	std::string s;
	std::string start = "[" + mon.getFactor().get_str() + "*" , end = "]", delim = "*";
    if (mon.getSize() > 0){
       	s += start;
       	for (int i = 0; i < mon.getSize(); i++){
       		if (this->phases[mon.vars[i]]) {
       			s += "x" + std::to_string(mon.vars[i]);
       		} else {
       			s += "f" + std::to_string(mon.vars[i]);
       		}
           	if (i != mon.getSize() - 1) s += delim;
       	}
       	s += end;
       	//Add factor at the end
    	}
   	else{
       	s += start;
       	s += end;
   	}
	return s;
}

//***************************************************************************************
std::string Polynom::to_string_opt() const {
    std::string start = "", end = "", delimPlus = "+", delimMinus = "-";
    std::string str;
    int size = this->getSet()->size();
    if (!this->getSet()->empty()){
        str += start;
        int num = 0;
        for (std::set<Monom>::iterator it=this->getSet()->begin(); it != this->getSet()->end(); ++it) {
			str += monToStringOpt(*it);
            num += 1;
            if (num != size) {
            	if ((std::next(it,1))->getFactor() >= 0) str += delimPlus;
            }
        }
    }
    else{
        str += "0";
    }
    return str;
}

//***************************************************************************************
std::string Polynom::to_string_with_phases_opt() const {
    std::string start = "", end = "", delimPlus = "+", delimMinus = "-";
    std::string str;
    int size = this->getSet()->size();
    if (!this->getSet()->empty()){
        str += start;
        int num = 0;
        for (std::set<Monom>::iterator it=this->getSet()->begin(); it != this->getSet()->end(); ++it) {
			str += monToStringWithPhasesOpt(*it);
            num += 1;
            if (num != size) {
            	if ((std::next(it,1))->getFactor() >= 0) str += delimPlus;
            }
        }
    }
    else{
        //s += start;
        //s += end;
        str += "0";
    }
    return str;
}

//***************************************************************************************
std::string Polynom::monToStringOpt(Monom mon) const {
	std::string s;
	std::string start = mon.getFactor().get_str() , delim = "*";
    if (mon.getSize() > 0){
       	s += start;
       	s += delim;
       	for (int i = 0; i < mon.getSize(); i++){
       		s += "x" + std::to_string(mon.vars[i]);
           	if (i != mon.getSize() - 1) s += delim;
       	}
    }
   	else{
       	s += start;
   	}
	return s;
}

//***************************************************************************************
std::string Polynom::monToStringWithPhasesOpt(Monom mon) const {
	std::string s;
	std::string start = mon.getFactor().get_str() , delim = "*";
    if (mon.getSize() > 0){
       	s += start;
       	s += delim;
       	for (int i = 0; i < mon.getSize(); i++){
       		if (this->phases[mon.vars[i]]) {
       			s += "x" + std::to_string(mon.vars[i]);
       		} else {
       			s += "f" + std::to_string(mon.vars[i]);
       		}
           	if (i != mon.getSize() - 1) s += delim;
       	}
    }
   	else{
       	s += start;
   	}
	return s;
}

//***************************************************************************************
void Polynom::replace_var_by_negation(std::string& str, varIndex var) {
	size_t pos = 0;
	std::string negString = "f";
	negString.append(std::to_string(var));
	std::string findString = "x";
	findString.append(std::to_string(var));
	str.replace(str.find(findString), findString.length(), negString);
}

//***************************************************************************************
std::vector<Monom> Polynom::modReductionWithQuotient(mpz_class modNum) {
	std::vector<Monom> toDelete;
	std::vector<Monom> quotient;
	for (auto& elem: this->polySet) {
		mpz_class coefBefore = elem.getFactor();
		mpz_mod(elem.factor.get_mpz_t(), elem.getFactor().get_mpz_t(), modNum.get_mpz_t());
		mpz_class coefAfter = elem.getFactor();
		if (coefBefore == coefAfter) continue;  // In this case no mod reduction was performed.
		if (coefAfter == 0) toDelete.push_back(elem);  // Mod reduced coef to 0. Remove this monomial.
		mpz_class diff = coefAfter - coefBefore;
		mpz_class fact = diff / modNum;
		quotient.push_back(elem);
		quotient.back().setFactor(fact);
	}
	for (auto& elem: toDelete) {
		this->eraseMonom(elem);
	}
	return quotient;
}

//***************************************************************************************
std::vector<std::string> Polynom::modReductionWithQuotientStr(mpz_class modNum) {
	std::vector<Monom> toDelete;
	std::vector<std::string> quotientStrVec;
	for (auto& elem: this->polySet) {
		mpz_class coefBefore = elem.getFactor();
		mpz_mod(elem.factor.get_mpz_t(), elem.getFactor().get_mpz_t(), modNum.get_mpz_t());
		mpz_class coefAfter = elem.getFactor();
		if (coefBefore == coefAfter) continue;  // In this case no mod reduction was performed.
		if (coefAfter == 0) toDelete.push_back(elem);  // Mod reduced coef to 0. Remove this monomial.
		mpz_class diff = coefAfter - coefBefore;
		mpz_class fact = diff / modNum;
		Monom tmpMon = elem;
		tmpMon.setFactor(fact);
		quotientStrVec.push_back(monToStringOpt(tmpMon));
	}
	for (auto& elem: toDelete) {
		this->eraseMonom(elem);
	}
	return quotientStrVec;
}

//***************************************************************************************
void Polynom::modReducePoly(mpz_class modNum) {
	std::vector<Monom> toDelete;
	for (auto& elem: this->polySet) {
		mpz_mod(elem.factor.get_mpz_t(), elem.getFactor().get_mpz_t(), modNum.get_mpz_t());
		if (elem.getFactor() == 0) toDelete.push_back(elem);
	}
	for (auto& elem: toDelete) {
		this->eraseMonom(elem);
	}
}

//***************************************************************************************
void Polynom::setModReduction(bool mode) {
	this->modReductionEnabled = mode;
}

//***************************************************************************************
void Polynom::setModReductionNumber(mpz_class modNum) {
	this->coefModReduction = modNum;
}

//***************************************************************************************
Monom Polynom::getShortestModel() {
	if (this->size() == 0) { std::cout << "Polynomial is empty. No model found." << std::endl; return Monom(); }
	int minSize = INT_MAX;
	Monom* monP = NULL;
	for (auto& elem: this->polySet) {
		if (elem.size < minSize)  {
			monP = const_cast<Monom*>(&elem);
			minSize = elem.size;
		}
	}
//	std::cout << "Shortest model of polynomial is: " << std::endl;
//	std::cout << *monP << std::endl;
	return *monP;
}

//***************************************************************************************
std::pair<std::string, std::string> Polynom::writeOutStartingPoly() {
	std::string startpoly = this->to_string();
	std::string modNumberStr = this->coefModReduction.get_str();
	return {startpoly, modNumberStr};
}

//***************************************************************************************
void Polynom::setProofGenerationMode(bool mode) {
	proofEnabled = mode; 	
}

//***************************************************************************************
void Polynom::startProofGeneration(std::string polyFile, std::string proofFile) {
	this->setProofGenerationMode(true);
	set_proof_filenames(polyFile, proofFile);
	writeStartPolyToFile(this->writeOutStartingPoly(), this->getVarSize());
}

//***************************************************************************************
std::string Polynom::writeReplacementAxiom(varIndex replace, std::list<Monom>& mons) {
	std::string returnStr = "-x";
	returnStr.append(std::to_string(replace));
	for (auto& elem: mons) {
		if (elem.getFactor() > 0) returnStr.append("+");
		returnStr.append(this->monToStringOpt(elem)); 
	}
	return returnStr;
}

//***************************************************************************************
std::string Polynom::writeReplacementAxiom(varIndex replace, std::set<Monom>* mons) {
	std::string returnStr = "-x";
	returnStr.append(std::to_string(replace));
	for (auto& elem: *mons) {
		if (elem.getFactor() > 0) returnStr.append("+");
		returnStr.append(this->monToStringOpt(elem)); 
	}

	return returnStr;
}

