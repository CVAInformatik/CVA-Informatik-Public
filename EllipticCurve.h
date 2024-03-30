#pragma once


#include "PrimeTable.h"
#include "CalculatorType.h"
#include "Calculator.h"
#include "Calculator2E30.h"


typedef std::vector<BINT*>  PrimeBaseListType;
typedef std::vector<BINT*>  InverseListType;

class  ResidueClass {

public:
	ResidueClass();
	~ResidueClass();

	void Setup( char* _Modulus, uint MaxFFTLength, uint ResidueLimit = 1000000);

	const PrimeBaseListType *PrimeBase;
	const InverseListType* Inverses;
	BINT* Modulus;

};

class PolynomialasResidue {

public:
	PolynomialasResidue(std::string& txt, ResidueClass res);
	~PolynomialasResidue();

	void AddTerm(std::string& txt, u64 exponent);
	void AddPolynomium(std::string& txt);


private:

	std::size_t MkTerms(std::size_t start, std::string text, std::string& constantExp, u64& exponent);
	std::size_t MkTermsConstant(std::size_t start, std::string text, std::string& constantExp, u64& exponent);
	std::size_t MkTermsExponent(std::size_t start, std::string text, std::string& constantExp, u64& exponent);

	void ParseFragment(std::vector<std::string>& fragments);

	std::vector<std::string> split(const std::string str, const std::string regex_str);

	ResidueClass Res;

	std::map<u64, BINT*>  PolyTerms;
};
