#pragma once

#include "CalculatorType.h"
#include "PrimeTable.h"
#include "Calculator.h"
#include "Calculator2E30.h"


//#ifdef CAL10
//#define CALCULATOR Calculator 
//#define BINT BInt
//#else
//#define CALCULATOR Calculator2E30 
//#define BINT BInt2E30
//#endif


void ModularExponentiation(BINT& Res, const BINT& a, const BINT& exp, const BINT& mod);

void ModularMultiplication(BINT& ab, const BINT& a, const BINT& b, const BINT& mod);

void ModularAddition(BINT& aplusb, const BINT& a, const BINT& b, const BINT& mod);

void ModularSquare(BINT& Res, const BINT& a, const BINT& mod);

void SquareRootModM(BINT& Res, BINT& A, BINT& M);

void SquareRootModPrime(BINT& Res, BINT& A, BINT& M);

bool MillerRabin(BINT& number, const std::vector<unsigned int>& witnesses);

void Factoring(char c[]);

void Faculty(BINT& res, int  a);

/* for experiments */

void Convert10E9to2E30(BInt2E30& dest, BInt& src);

void Convert2E30to10E9(BInt& dest, BInt2E30& src);

void MersenneBInt2E20(BInt2E30& dest, uint N);


std::string* ItoA(Calculator& cal);

std::string* ItoA(Calculator2E30& cal);
