#pragma once


/*
Copyright  © 2024 Claus Vind - Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/


//
//CALCULATOR::CALCULATOR() {
//	dist = new std::uniform_int_distribution<uint>(0, RMOD - 1);
//}
//
//CALCULATOR::~CALCULATOR() {
//	delete dist;
//}

void CALCULATOR::Dup() { 
	if (stack.size() > 0) {
		BINTPTR temp = stack.back();
		stack.push_back(temp);
	}
};  // push a copy of TOS on Stack


bool CALCULATOR::IsEven()
{
	if (stack.size() > 0) 	return ((stack.back()->number[0] & 1) == 0);
	return false;
}

bool CALCULATOR::IsEqual(int n)
{
	if (stack.size() > 0) 	
		return  (stack.back()->number.size() ==1 ) && (n == (stack.back()->number[0] * (stack.back()->sign)));
	return false;
}

/* does not require a normalized BINT ! */
bool CALCULATOR::IsZero(const BINT &arg)
{
	int num = 0;
	for (unsigned int i = 0; !num && (i < arg.number.size()); i++)	num |= arg.number[i];
	return num == 0;
}


bool CALCULATOR::IsEqual() {
	size_t sz = stack.size();
	if (sz < 2)
	{
		std::cout << "IsEqual(): not enough arguments" << std::endl;
		return false;
	}
	else
		return IsEqual(*stack[sz - 1], *stack[sz - 2]);
	return true;
}

bool CALCULATOR::IsEqual(BINT A, BINT B)
{
	if (A.sign != B.sign) return false;
	if (A.number.size() != B.number.size()) return false;
	for (size_t i = 0; i < A.number.size(); i++)
		if (A.number[i] != B.number[i]) return false;
	return true;
}


bool CALCULATOR::IsAbiggerNummerically(BINTPTR A, BINTPTR B) 
{ 
	return IsAbiggerNummerically(*A, *B); 
}


uint  CALCULATOR::_Rand(uint UpperBound)
{
	uint ix2 = dist->operator()(rd);
	return ix2 % UpperBound;
}



void CALCULATOR::Swap()
{

	if (stack.size() > 1) {

		BINTPTR t1 = stack.back(); 		stack.pop_back();
		BINTPTR t2 = stack.back();		stack.pop_back();
		stack.push_back(t1);
		stack.push_back(t2);
	}
}

void CALCULATOR::Push(char* c) {
	BINT temp;
	char* ct1, * ct2;

	ct1 = ct2 = c;
	temp.sign = 1;

	while (isspace(*ct1)) ct1++;
	if (*ct1 == '-') {
		temp.sign = -1;
		ct1++;
	}
	ct2 = ct1;

	while (isdigit(*ct2)) ct2++; // find the end;
	if (ct2 == ct1) {  //bail out something is wrong
		std::cerr << "unknown format " << c << std::endl;
	}
	// we assume we have something number like 
	temp.number.clear();
	temp.number.push_back(0);
	int tempInt = 0;
	int counter = 0;
	while (isdigit(*ct1)) {
		tempInt = 10 * tempInt;
		tempInt += (*ct1) - '0';
		counter++;
		ct1++;
		if (counter == 3) {
			Push(tempInt);
			SimpleAdditionSubtractionLadder1(1, 1000, 1, temp);
			Add();
			Pop(temp);
			counter = tempInt = 0;
		}
	}
	switch (counter) {
	case 1:
		Push(tempInt);
		SimpleAdditionSubtractionLadder1(1, 10, 1, temp);
		Add();
		Pop(temp);
		break;
	case 2:
		Push(tempInt);
		SimpleAdditionSubtractionLadder1(1, 100, 1, temp);
		Add();
		Pop(temp);
		break;
	default: break;
	}
	Normalize(temp);
	Push(temp);
}


void CALCULATOR::Push(int i) {
	BINTPTR temp(new BINT); temp->number.clear();
	s64 it = i;
	temp->sign = (it >= 0) ? 1 : -1;

	if (it < 0) it = -1 * i;
	do {
		temp->number.push_back(it % (s64)RMOD);
		it = it / RMOD;
	} while (it > 0);
	stack.push_back(temp);
}


void CALCULATOR::Push(const BINT& b)
{
	BINTPTR temp(new BINT);
	Dup(*temp, b);
	stack.push_back(temp);
}

void CALCULATOR::ChangeSign() {
	if (stack.size() > 0) {
		BINTPTR temp(new BINT); Pop(*temp);
		temp->sign = -1 * temp->sign;
		stack.push_back(temp);
	}

}


int CALCULATOR::TOSSize() {
	if (stack.size() > 0) return (int)stack.back()->number.size();
	return 0;
}


/* other */
int  CALCULATOR::StackSize() { return  (int)stack.size(); }; // how deep is the stack...
int  CALCULATOR::StoreSize() { return  (int)Store.size(); }; // how many items in the store...


/* store operations  */

void CALCULATOR::ClearStore(const std::string& loc) {
	if (Store.find(loc) != Store.end())
		Store.erase(loc);
	else
		std::cout << "Lookup failed" << std::endl;
}


void CALCULATOR::PushStore(const std::string& loc)
{
	if (Store.find(loc) != Store.end())
		stack.push_back(Store.find(loc)->second);
	else
		std::cout << "Lookup failed" << std::endl;
}


void CALCULATOR::PopStore(const std::string& loc)
{

	if (stack.size() > 0) {
		Store[loc] = stack.back();
		stack.pop_back();
	}
};

void CALCULATOR::Dup(BINT& D, const BINT& S)
{
	D.number.clear();
	D.sign = S.sign;
	for (int ix = 0; ix < S.number.size(); ix++)
		D.number.push_back(S.number[ix]);
}

void CALCULATOR::Normalize(BINT& b) {
	while (b.number.size() && b.number.back() == 0) b.number.pop_back();
	if (b.number.size() == 0) {
		b.number.push_back(0);
		b.sign = 1;
	}
}


void CALCULATOR::Mul2(BINT& A) {
	int carry = 0;

	for (int ix = 0; ix < A.number.size(); ix++)
	{
		int t = (A.number[ix] << 1) + carry;
		carry = 0;
		if (t >= RMOD) {
			carry = 1;
			t -= RMOD;
		}
		A.number[ix] = t;
	}

	if (carry) A.number.push_back(carry);
}


void CALCULATOR::Add() {
	if (stack.size() < 2)
		std::cout << "Add() needs two arguments" << std::endl;
	else {
		BINTPTR temp(new BINT);

		Dup(*temp, *stack.back());
		stack.pop_back();
		AddAux(temp->sign, *temp, stack.back()->sign, *stack.back());
		stack.pop_back();
		Normalize(temp);
		stack.push_back(temp);
	}
}



void CALCULATOR::SimpleAdditionSubtractionLadder1(int sign, s64 A, int BSign, const BINT& B)
{
	{
		/* this is taken from Crandall& Pomerance
		*   "Prime Numbers,  A Computational Perspective" 2nd edition
		*/
		BINTPTR result(new BINT); result->number.clear();result->number.push_back(0);

		s64 A3 = A * 3;
		s64 A1 = A;
		s64 A3EXORA1 = A3 ^ A1;
		s64 Mask = (A3 >> 1 | A3 >> 2); // ! 
		Mask |= Mask >> 2;
		Mask |= Mask >> 4;
		Mask |= Mask >> 8;
		Mask |= Mask >> 16;
		Mask |= Mask >> 32;
		Mask++; // is now one '1' positioned the first non-zero bit of A3

		Dup(*result, B);
		Mask = Mask >> 1;
		while (Mask > 1) {
			Mul2(*result);
			if (Mask & (A3EXORA1)) {
				if (Mask & A3) {
					AddAux(1, *result, 1, B);
				}
				else {
					AddAux(1, *result, -1, B);
				}
			}
			Mask = Mask >> 1;
		}
		Normalize(*result);
		if (!IsZero(*result)) result->sign = sign;

		stack.push_back(result);
	}

}



void CALCULATOR::Div2(BINT& A) {

	int borrow = 0;
	for (s64 ix = A.number.size() - 1; ix >= 0; ix--)
	{
		int t = A.number[ix] + borrow; borrow = 0;
		if (t & 1)  borrow = RMOD;
		A.number[ix] = t >> 1;
	}
	while (A.number.size() && A.number.back() == 0) A.number.pop_back();
}

void CALCULATOR::Exp()
{
	if (stack.size() < 2)
		std::cout << "Exp() needs two arguments" << std::endl;
	else {
		BINT    Exponent;  Pop(Exponent);
		BINTPTR Argument(new BINT); Pop(*Argument);
		BINTPTR  Result(new BINT); Result->number.push_back(1);
		stack.push_back(Result);
		while (!IsZero(Exponent)) {
			if (Exponent.number.size() && Exponent.number[0] & 1) {
				stack.push_back(Argument);
				Mul();
			}
			Div2(Exponent);
			if (!IsZero(Exponent)) {
				stack.push_back(Argument);
				stack.push_back(Argument);
				Mul();
				Dup(*Argument, *stack.back());
				stack.pop_back();
			}
		}
	}
}

#define SIMPLEMULT(s,A,sb, B)  SimpleAdditionSubtractionLadder1(s,A,sb, B)

void CALCULATOR::RussianPeasantMult() {
	BINT x;
	BINT y;

	Pop(x); 	Pop(y);
#if TESTMUL == 1
	std::cout << "RussianPeasantMult() " << std::endl;
	DumpInt(" x: ", x);
	DumpInt(" y: ", y);
#endif
	if (IsZero(x) || IsZero(y))
	{
		BINTPTR result(new BINT());
		result->number.push_back(0);
		result->sign = 1;
		stack.push_back(result);

	}
	else {
		s64 xint = 0;
		s64 yint = 0;

		if (x.number.size() < 1 + SMALLNUMBERLIMIT)
			for (s64 ix = x.number.size(); ix > 0;ix--)
				xint = (xint * RMOD) + x.number[ix - 1];
		if (y.number.size() < 1 + SMALLNUMBERLIMIT)
			for (s64 iy = y.number.size(); iy > 0;iy--)
				yint = (yint * RMOD) + y.number[iy - 1];

		if (yint == 0)  		SIMPLEMULT(x.sign * y.sign, xint, y.sign, y);
		else if (xint == 0)     SIMPLEMULT(y.sign * x.sign, yint, x.sign, x);
		else if (xint < yint) 	SIMPLEMULT(x.sign * y.sign, xint, y.sign, y);
		else                	SIMPLEMULT(y.sign * x.sign, yint, x.sign, x);
	}
}

/*

the code below is adapted from

		Algorithm 14.9 	in

		The Handbook of Applied Cryptograpy by A.Menezes, P van Oorschot and S.Vanstone (CRC Press 1996).
*/

/*
  The Asign and Bsign parameters are normally A and B real signs, but sometimes we
  want to overrule that.
*/

void CALCULATOR::AddAux(int Asign, BINT& A, int Bsign, const BINT& B)
{
	/* setup  */
	u64 min_sz = std::min(A.number.size(), B.number.size());

	if (Asign == Bsign)
	{
		/* they have identical signs */
		/* this is done in twp steps
		*    1)  Add the digits without carry
		*    2)  do ' + carry mod RMOD'
		*/

		/* step 1: add the other number, it may be longer or shorter than temp */
		for (u64 i = 0; i < min_sz; i++)  A.number[i] = A.number[i] + B.number[i];

		for (u64 i = min_sz; i < B.number.size(); i++)	A.number.push_back(B.number[i]);

		/* step 3 */
		int carry = 0;
		for (u64 i = 0; i < A.number.size(); i++)
		{
			A.number[i] += carry;
			if (A.number[i] >= RMOD) {
				A.number[i] -= RMOD;
				carry = 1;
			}
			else
				carry = 0;
		}
		if (carry) A.number.push_back(carry);
	}
	else
	{
		/* they have different signs */
		/* this is done in four steps
		*    1)  Multiply the digits in the negative number with -1
		*    2)  Add the digits without carry
		*    3)  do ' + carry mod RMOD'
		*    if carry out from the MSD is -1
		*    4)   (optional)  subtract from 0 and adjust sign of result.
		*/

		/* step 1: negate digits in the negative number, which is TOS,  and save them in temp  */
		if (Asign == -1) {
			A.sign = 1;
			for (u64 i = 0; i < A.number.size(); i++) A.number[i] = -1 * A.number[i];
			for (u64 i = 0; i < min_sz; i++)  A.number[i] = A.number[i] + B.number[i];
			for (u64 i = min_sz; i < B.number.size(); i++)	A.number.push_back(B.number[i]);
		}
		else { /* Bsign == -1 */
			for (u64 i = 0; i < min_sz; i++)  A.number[i] = A.number[i] - B.number[i];
			for (u64 i = min_sz; i < B.number.size(); i++)	A.number.push_back(-B.number[i]);
		}

		/* step 3: now adjust carry/borrow and mod */
		int c = 0;
		for (int i = 0; i < A.number.size(); i++)
		{
			A.number[i] += c;
			if (A.number[i] < 0) {
				c = -1;
				A.number[i] += RMOD;
			}
			else
				c = 0;
		}

		/* step 4: final correction if c is set */
		if (c != 0) {
			A.sign = -1;
			int c1 = 0;
			for (int i = 0; i < A.number.size(); i++)
			{
				A.number[i] = 0 - A.number[i] + c1;
				if (A.number[i] < 0) {
					c1 = -1;
					A.number[i] += RMOD;
				}
				else
					c1 = 0;
			}
		}
	}
}


void CALCULATOR::Jacobi()
{
	if (stack.size() < 2)
		std::cout << "Jacoby/Legendre  needs two arguments" << std::endl;
	else {
		BINTPTR Res(new BINT); Res->number.clear();
		BINTPTR A(new BINT); Pop(*A);
		BINTPTR M(new BINT); Pop(*M);
		if ((M->number.size() == 1) && (M->number[0] == 1))
		{
			Res->number.clear(); Res->number.push_back(1); Res->sign = 1;
			stack.push_back(Res);
			return;
		}
		else {

			stack.push_back(A);
			stack.push_back(M);
			Mod();
			if (IsZero(*stack.back())) {  //A|M 
				stack.pop_back();
				Res->number.push_back(0); Res->sign = 1;
				stack.push_back(Res);
			}
			else {
				Pop(*A); //reminder of A mod M
				Res->number.push_back(1); Res->sign = 1;
				while (!IsZero(*A)) {
					while ((A->number[0] & 1) == 0)
					{
						Div2(*A);
						switch (M->number[0] & 0x7)
						{
						case 3: case 5:   Res->sign = -1 * Res->sign;
							break;
						default:
							break;
						}
					}
					BINT temp;
					Dup(temp, *A);
					Dup(*A, *M);
					Dup(*M, temp);
					if ((3 == (A->number[0] & 0x3)) && (3 == (M->number[0] & 0x3)))
						Res->sign = -1 * Res->sign;
					stack.push_back(A);
					stack.push_back(M);
					Mod();
					Pop(*A);
				}
				if ((M->number.size() == 1) && (M->number[0] == 1))  stack.push_back(Res);
				else
				{
					Res->number.clear(); Res->number.push_back(0); Res->sign = 1;
					stack.push_back(Res);
					return;
				}
			}
		}
	}
}




void CALCULATOR::GCDAux(BINTPTR X, BINTPTR Y)
{
	BINT g;
	BINT x;  		Dup(x, *X);
	BINT y;  		Dup(y, *Y);

	/* remove common 1000 factors */
	u64 min_sz = std::min(x.number.size(), y.number.size());
	int i = 0;
	for (u64 i = 0; i < min_sz; i++)
		if ((x.number[i] | y.number[i]) == 0)  i++; /* counter #1000 */
		else break;

	/* then remove them and update g */
	for (int i2 = 0; i2 < i; i2++)  /* then remove them and update g */
	{
		g.number.push_back(0);
		DIVRMOD(x); DIVRMOD(y);
	}
	g.number.push_back(1);


	/* copy to temps*/
	BINT u; Dup(u, x);
	BINT v; Dup(v, y);

	/* remove common 2 factors */
	while ((u.number[0] & 1) == 0 && (v.number[0] & 1) == 0)
	{
		Div2(u); Div2(v); Mul2(g);
	};

	BINT A; A.number.push_back(1);
	BINT B; B.number.push_back(0);
	BINT C; C.number.push_back(0);
	BINT D; D.number.push_back(1);

	while (!IsZero(u)) {
		//std::cout << "-----------------------------------" << std::endl;
		DUMPINT("u", u); 		DUMPINT("v", v);		DUMPINT("A", A);
		DUMPINT("B", B);		DUMPINT("C", C);		DUMPINT("D", D);
		//std::cout << "-----------------------------------" << std::endl;

		while ((u.number[0] & 1) == 0)
		{
			Div2(u); DUMPINT("u", u);
			if ((A.number[0] & 1) == 0 && (B.number[0] & 1) == 0) {
				DUMPINT("A", A);
				Div2(A); DUMPINT("B", B);
				Div2(B);
			}
			else {
				DUMPINT("A", A);
				Add(A, y);		DUMPINT("A", A);
				Div2(A);		DUMPINT("A", A);DUMPINT("B", B);
				Sub(B, x);		DUMPINT("B", B);
				Div2(B);		DUMPINT("B", B);
			}
		}

		while ((v.number[0] & 1) == 0) {
			DUMPINT("v", v);
			Div2(v); DUMPINT("v", v);
			if ((C.number[0] & 1) == 0 && (D.number[0] & 1) == 0) {
				DUMPINT("C", C);
				Div2(C);  DUMPINT("D", D);
				Div2(D);
			}
			else {
				DUMPINT("C", C);  DUMPINT("y", y);
				Add(C, y);	DUMPINT("C", C);
				Div2(C);	DUMPINT("C", C);  DUMPINT("D", D);	DUMPINT("x", x);
				Sub(D, x);	DUMPINT("D", D);
				Div2(D);	DUMPINT("D", D);
			}
		}
		DUMPINT("u", u);	DUMPINT("v", v);

		if (IsALarger(u, v) || IsEqual(u, v))
		{
			DUMPINT("u", u);DUMPINT("v", v);
			Sub(u, v);	DUMPINT("u", u);	DUMPINT("A", A);	DUMPINT("C", C);
			Sub(A, C);	DUMPINT("A", A);	DUMPINT("B", B);	DUMPINT("D", D);
			Sub(B, D);	DUMPINT("B", B);
		}
		else {
			DUMPINT("u", u);   DUMPINT("v", v);
			Sub(v, u);	DUMPINT("v", v);  DUMPINT("A", A);	DUMPINT("C", C);
			Sub(C, A);	DUMPINT("C", C);  DUMPINT("B", B);	DUMPINT("D", D);
			Sub(D, B);	DUMPINT("D", D);

		}
	};


	BINTPTR a1(new BINT); Dup(*a1, C);
	BINTPTR b1(new BINT); Dup(*b1, D);
	BINTPTR gcd(new BINT); Dup(*gcd, v);
	BINTPTR gmul(new BINT); Dup(*gmul, g);

	stack.push_back(a1);
	stack.push_back(b1);
	stack.push_back(gcd);
	stack.push_back(gmul);
	Mul();
}

void CALCULATOR::GCD() {
	if (stack.size() < 2)
		std::cout << "GCD() needs two arguments" << std::endl;
	else {
		if (IsZero(*stack[stack.size() - 1]))
			std::cout << "GCD() needs non-zero arguments" << std::endl;
		else if (IsZero(*stack[stack.size() - 2]))
			std::cout << "GCD() needs non-zero arguments" << std::endl;
		else {
			BINTPTR A = stack.back(); stack.pop_back();
			BINTPTR B = stack.back(); stack.pop_back();
			GCDAux(A, B);
		}
	}


}

void CALCULATOR::dumpStack(int p)
{
	printf("stack(%d): \n", p);
	if (stack.size() == 0)
		std::cout << "empty stack" << std::endl;
	else
		for (int ix = 1; ix <= stack.size();ix++) {
			printf("level %-2d : ", ix - 1);
			DUMPINT("", *stack[stack.size() - ix]);
		}
	printf("--------------------------\n");
}