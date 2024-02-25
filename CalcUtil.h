#pragma once

#include "PrimeTable.h"
#include "Calculator.h"

void ModularExponentiation(BInt& Res, const BInt &a, const BInt &exp, const BInt &mod);

void ModularMultiplication(BInt& ab, const BInt &a, const BInt &b, const BInt &mod);

void ModularAddition(BInt& aplusb, const BInt &a, const BInt &b, const BInt &mod);

void ModularSquare(BInt& Res, const BInt &a, const BInt &mod);

void SquareRootModPrime(BInt& Res, BInt& A, BInt& M);

void SquareRootModM(BInt& Res, BInt& A, BInt& M);

bool MillerRabin(BInt& number, const std::vector<unsigned int>& witnesses);

void Factoring(char c[]);

void Faculty(BInt& res, int  a);

/* for experiments */
typedef  struct _BInt2E30
{
	int sign;
	std::vector<int>  number;
	_BInt2E30() { sign = 1; number.clear(); };
} BInt2E30;


void Convert10E9to2E30(BInt2E30& dest, BInt& src);

void Convert2E30to10E9(BInt& dest, BInt2E30& src);

void MersenneBInt2E20(BInt2E30& dest, uint N);
