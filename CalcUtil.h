/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/

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
