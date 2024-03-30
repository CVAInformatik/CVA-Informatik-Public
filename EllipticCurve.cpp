
#include "PrimeTable.h"

#include "EllipticCurve.h"

#include <regex>

ResidueClass::ResidueClass()
{
	PrimeBase = NULL;
	Modulus = NULL;
}

ResidueClass::~ResidueClass()
{
	if (PrimeBase != NULL) delete PrimeBase;
	if (Modulus != NULL) delete Modulus;
}

void ResidueClass::Setup(char* _Modulus, uint MaxFFTLength, uint ResidueLimit)
{
	PrimeTable pt(ResidueLimit);

	CALCULATOR cal;

	//Modulus = new BINT();
	PrimeBaseListType *tPrimeBase = new PrimeBaseListType();
	InverseListType *tInverses = new InverseListType();

	cal.Push( _Modulus);
	cal.PopStore("Modulus");
	cal.Push(1);
	cal.PopStore("Product");

	uint Candidate = ResidueLimit;
	if (0 ==  (Candidate % 2)) Candidate--;
	for (; Candidate > 2; Candidate -= 2) {
		cal.Clear(); //we start we a clean stack
		if (!pt.IsPrime(Candidate)) continue;
		cal.PushStore("Modulus");
		cal.Push(Candidate);
		cal.Dup();
		BINT* cand = new BINT();
		cal.Pop(*cand);
		cal.PushStore("Product");
		cal.Push(*cand);
		tPrimeBase->push_back(cand);
		cal.Mul();
		cal.PopStore("Product");
		cal.PushStore("Modulus");
		cal.PushStore("Product");
		if (cal.IsLarger()) break;
		/* now we have a primebase which multiplied to gether is larger than modulus */
	}

	Inverses = tInverses;
	PrimeBase = tPrimeBase;
}

/* 
		given a string like " nnnnnX^E1 + nnnnnX^E2+ nnnnnX^E3....."
		and a ResidueClass instance, a CRT represenation of the polynomium is created.

*/
PolynomialasResidue::PolynomialasResidue(std::string& txt, ResidueClass res) 
{


	Res = res;

	AddPolynomium(txt);
}



/* terms  of a polynomial and returns the position it stopped, for the next term

	The formal variable is always (case-independent)  'X'/'x'
	Exponents are always positive !

 The following output for various terms

		input					constantExp    exponent	
        "x +....   "			"1",				1
		" -62834682 +......		"-6262834682"		0
		" +846283468X^25 +.....	"+846283468"		25
		syntax error in string  ""                  0	

*/

std::size_t PolynomialasResidue::MkTermsExponent(std::size_t start, std::string text, std::string& constantExp, u64& exponent)
{
	exponent = 0;
	size_t pos = start;

	for (;pos < text.size(); )
		switch (text[pos])
		{
		case '0': case '1':	case '2': case '3':	case '4':  case '5': case '6': 	case '7': case '8':	case '9':
			exponent = exponent * 10;
			exponent = exponent + text[pos] -'0';
			pos++;
			break;
		case '+': 	case '-':	case ' ':	case '\t':	return pos;
			break;
		default:
			constantExp = "";
			exponent = 0;
			return pos;
			break;
		}

	constantExp = "";
	exponent = 0;
	return pos;
}

std::size_t PolynomialasResidue::MkTermsConstant(std::size_t start, std::string text, std::string& constantExp, u64& exponent)
{
	size_t pos = start;

	constantExp = "";
	for(;pos< text.size(); )
	switch (text[pos])
	{
	case '+': case '-': case '0': case '1':	case '2': case '3': case '4': case '5': case '6': case '7': case '8': 	case '9':
		constantExp = constantExp.append(1, text[pos]); pos++;
		break; 
	case ' ': 	case '\t': pos++;
		break;
	case 'x': 	case 'X': 
		if (pos + 1 < text.size() && text[pos + 1] == '^')
			return MkTermsExponent(pos + 2, text, constantExp, exponent);
		else {
			exponent = 1;
			return pos + 1;
		}
		break;
	default:
		constantExp = "";
		exponent = 0;
		return pos;
		break;
	}
	constantExp = "";
	exponent = 0;
	return pos;
}


std::size_t PolynomialasResidue::MkTerms(std::size_t start, std::string text, std::string& constantExp, u64& exponent)
{
	size_t pos = start;

	while (pos < text.size() && isspace(text[pos])) pos++;
	if (pos < text.size()) return MkTermsConstant(pos, text, constantExp, exponent);
	constantExp = "";
	exponent = 0;
	return  pos;
}


void PolynomialasResidue::AddTerm(std::string& txt, u64 exponent)
{
	CALCULATOR   Cal;
	if (PolyTerms.find(exponent) == PolyTerms.end())
	{  /* new exponent, we add a 0 entry for the constant */
		BINT* b = new BINT();
		Cal.Push(0);
		Cal.Pop(*b);
		PolyTerms[exponent] = b;
	}
	Cal.Push((char*)txt.c_str());
	Cal.Push(*PolyTerms[exponent]);
	Cal.Add();
	BINT* b1 = new BINT();
	Cal.Pop(*b1);
	PolyTerms[exponent] = b1;

}

void PolynomialasResidue::AddPolynomium(std::string& txt)
{

	std::string Const;
	u64          Exp;
	CALCULATOR   Cal;

	size_t pos = 0;
	while (pos < txt.size())
	{
		pos = MkTerms(pos, txt, Const, Exp);

		if (Const == "" && Exp == 0)  return;

		if (PolyTerms.find(Exp) == PolyTerms.end())
		{  /* new exponent, we add a 0 entry for the constant */
			BINT* b = new BINT();
			Cal.Push(0);
			Cal.Pop(*b);
			PolyTerms[Exp] = b;
		}
		Cal.Push((char *) Const.c_str());
		Cal.Push(*PolyTerms[Exp]);
		Cal.Add();
		BINT* b1 = new BINT();
		Cal.Pop(*b1);
		PolyTerms[Exp] = b1;
	}

}



PolynomialasResidue::~PolynomialasResidue(){
}

