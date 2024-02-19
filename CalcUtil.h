#pragma once

#include "PrimeTable.h"
#include "Calculator.h"

void ModularExponentiation(BInt& Res, const BInt &a, const BInt &exp, const BInt &mod);

void ModularMultiplication(BInt& ab, const BInt &a, const BInt &b, const BInt &mod);

void ModularAddition(BInt& aplusb, const BInt &a, const BInt &b, const BInt &mod);

void ModularSquare(BInt& Res, const BInt &a, const BInt &mod);

void SquareRootModPrime(BInt& A, BInt& M);

void SquareRootModM(BInt& A, BInt& M);

bool MillerRabin(BInt& number, const std::vector<unsigned int>& witnesses);