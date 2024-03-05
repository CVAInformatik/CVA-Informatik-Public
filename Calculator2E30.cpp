/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/


#include "Calculator2E30.h"



Calculator2E30::Calculator2E30() {
	dist = new std::uniform_int_distribution<uint>(0, RMOD - 1);

}

Calculator2E30::~Calculator2E30() {
	delete dist;
}
// stack operations
void Calculator2E30::Push(int i) {
	BInt2E30Ptr temp(new BInt2E30); temp->number.clear();
	s64 it = i;
	temp->sign = (it >= 0) ? 1 : -1;

	if (it < 0) it = -1 * i;
	do {
		temp->number.push_back(it % (s64)RMOD);
		it = it / RMOD;
	} while (it > 0);
	stack.push_back(temp);

}
void Calculator2E30::Push(const BInt2E30& b) {
	BInt2E30Ptr temp(new BInt2E30);
	Dup(*temp, b);
	stack.push_back(temp);

}
void Calculator2E30::Swap() {

	if (stack.size() > 1) {
		BInt2E30Ptr t1 = stack.back(); 		stack.pop_back();
		BInt2E30Ptr t2 = stack.back();		stack.pop_back();
		stack.push_back(t1);
		stack.push_back(t2);
	}
}

/* ASCII Conversions    */
void Calculator2E30::Push(char* c) {
	BInt2E30 temp;
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

/* Arithmetic */

void Calculator2E30::Mul() {
#define OVERALLOCATION 2
	if (stack.size() < 2)
		std::cout << "Mul() needs two arguments" << std::endl;
	else if (SMALLNUMBERLIMIT && stack[stack.size() - 1]->number.size() < 1 + SMALLNUMBERLIMIT)
		RussianPeasantMult();
	else if (SMALLNUMBERLIMIT && stack[stack.size() - 2]->number.size() < 1 + SMALLNUMBERLIMIT)
		RussianPeasantMult();
	else {
		PrimeFactorDFT pf;
		BInt2E30Ptr temp(new BInt2E30); temp->number.clear();

		temp->sign = stack[stack.size() - 1]->sign * stack[stack.size() - 2]->sign;

		factorSeq  factors;

		u64 min_sz = stack[stack.size() - 1]->number.size() +
			stack[stack.size() - 2]->number.size();

		u64 Length = 1;

		int factorlist[] = { 31,19,17,13,11,7,5,3,2 };

		for (int i = 0; i < sizeof(factorlist) / sizeof(factorlist[0]); i++)
		{
			factors.push_back(factorlist[i]);
			Length *= factors.back();
			if (Length > 3 * 2 * min_sz) break; /* 3 because we go from Radix 10^9 to Radix 10^3 */
			if (i == (sizeof(factorlist) / sizeof(factorlist[0]) - 1)) {
				std::cout << "Sorry Number too Large" << std::endl;
				return;
			}
		}

		pf.SetFactors(factors);

		if (pf.Status() > 0) {
			/* this should be less dynamic, but right now it is OK */
			Data* real1 = new Data[pf.Status() + OVERALLOCATION];
			Data* imag1 = new Data[pf.Status() + OVERALLOCATION];
			Data* real3 = new Data[pf.Status() + OVERALLOCATION];
			Data* imag3 = new Data[pf.Status() + OVERALLOCATION];
			for (s64 i = 0; i < pf.Status() + OVERALLOCATION; i++) {
				real1[i] = 0;
				imag1[i] = 0;
				real3[i] = 0;
				imag3[i] = 0;
			}

			LoadFFT(stack.back(), real1);
			stack.pop_back();
			LoadFFT(stack.back(), imag1);
			stack.pop_back();

			pf.forwardFFT(real1, imag1);

			real3[0] = real1[0] * imag1[0];
			s64 size = pf.Status();
			for (s64 i = 1; i < size; i++)
			{
				Data X01Real = real1[i];
				Data X01Imag = imag1[i];
				Data X02Real = real1[size - i];
				Data X02Imag = imag1[size - i];
				Data X1Real = (X01Real + X02Real) / 2.0;
				Data X1Imag = (X01Imag - X02Imag) / 2.0;
				Data X2Imag = -1 * (X01Real - X02Real) / 2;
				Data X2Real = (X01Imag + X02Imag) / 2;
				Data X3Real = X1Real * X2Real - X1Imag * X2Imag;
				Data X3Imag = X1Real * X2Imag + X1Imag * X2Real;
				real3[i] = X3Real;
				imag3[i] = X3Imag;
			}
			for (s64 i = 0; i < size; i++)
			{
				real1[i] = real3[i];
				imag1[i] = imag3[i];
			}

			pf.ScaledInverseFFT(real1, imag1);

			Carry(size, real1);

			/* convert back to radix 2^30 from radix 2^10 double*/

			for (int i = 0; i < pf.Status(); i += 3)
			{
				int t0 = (int)real1[i];
				int t1 =  1024 * (int)real1[i + 1];  // these values are 0 when we reach
				int t2 = 1024 * 1024 * (int)real1[i + 2];// the end of buffer, due to
				temp->number.push_back(t0 + t1 + t2);//OVERALLOCATION
			}

			Normalize(temp);

			stack.push_back(temp);

			delete[] real1;
			delete[] imag1;
			delete[] real3;
			delete[] imag3;
		}
	}
}

#define RMOD2E30 0x3FF

void Calculator2E30::LoadFFT(BInt2E30Ptr A, Data* Buffer)
{
	int carry = 0;
	int FFTIndex = 0;
	for (int ix = 0; ix < A->number.size(); ix++) {
		int temp = A->number[ix];   // temp is  radix 2^30
		for (int i = 0; i < 3; i++) {
			int tmp =  temp & RMOD2E30;// tmp is  radix 10^3
			temp = temp >> 10 ;
			tmp = tmp + carry; 
			carry = 0;
			if (tmp >= (RMOD2E30+1)/2 ) {  // tmp is 'balanced' radix 10^3 
				tmp = tmp - (RMOD2E30+1);
				carry++;
			}
			Buffer[FFTIndex++] = (double)tmp;
		}
	}
	if (carry)
		Buffer[FFTIndex] = 1.0;
}

void Calculator2E30::Carry(s64 size, Data* Buffer)
{
	s64 carry = 0;
	s64 tmp = 0;

	/* conversion from 'balanced' Radx 2^10  double  to unbalanced Radix 2^10 double */
	for (s64 i = 0; i < size; i++)
	{
		tmp = (s64)std::round(Buffer[i]);

		tmp += carry; carry = 0;
		while (tmp < (-1 * ((RMOD2E30 +1)/ 2))) { tmp += (RMOD2E30+1); carry--; }
		while (tmp >= (((RMOD2E30 +1)/ 2) )) { tmp -= (RMOD2E30+1); carry++; }
		Buffer[i] = (double)tmp;
	}
	carry = 0;
	/* carry we are still in Radix 2^10 double*/
	for (s64 i = 0; i < size; i++)
	{
		tmp = (s64) std::round(Buffer[i]);
		tmp += carry;  carry = 0;
		if (tmp < 0) { tmp += (RMOD2E30+1); carry--; }
		Buffer[i] = double(tmp);
	}

}


void Calculator2E30::QuotientRemainder() {
	BInt2E30Ptr Quotient(new BInt2E30); Quotient->number.clear();
	Quotient->sign = 1;

	int counter = 0;


	if (stack.size() < 2)
		std::cout << "QuotientReminder() needs two arguments" << std::endl;
	else if (IsZero(*stack[stack.size() - 1])) {
		std::cout << "divison by zero" << std::endl;
		return;
	}
	else if ((stack[stack.size() - 1]->number.size() == 1) &&
		(stack[stack.size() - 2]->number.size() == 1)) {
		/* small numbers both less than RMOD */
		BInt2E30Ptr    Remainder(new BInt2E30); Remainder->number.clear();
		Quotient->number.push_back(stack[stack.size() - 2]->number[0] / stack[stack.size() - 1]->number[0]);
		Remainder->number.push_back(stack[stack.size() - 2]->number[0] % stack[stack.size() - 1]->number[0]);
		stack.pop_back(); stack.pop_back();
		stack.push_back(Quotient);
		stack.push_back(Remainder);
	}
	else {

		//dumpStack(0);
		BInt2E30    divisor;  Pop(divisor);
		BInt2E30    dividend; Pop(dividend);
		DUMPINT("divisor  ", divisor);
		DUMPINT("dividend ", dividend);

		/* simple first draft !*/
			/* pseudo code
			Make an approximation to 1/divisor  : reciprocal
			Quotient = dividend *reciprocal
			do
				form dividend -  Quotient * divisor  =  reminder  // if negative Quotient too big ;
				if(reminder > 0 and reminder < divisor )
					we are done.
				else
					Quotient = Quotient + reminder *reciprocal ;
			while  true
		*/
		int     reciprocal = RMOD / (2+divisor.number.back());
		BInt2E30Ptr _reciprocal(new BInt2E30); _reciprocal->number.clear(); _reciprocal->number.push_back(reciprocal);
		int     shift = (int)divisor.number.size();
		BInt2E30Ptr _dividend(new BInt2E30); Dup(*_dividend, dividend);
		BInt2E30Ptr _divisor(new BInt2E30);  Dup(*_divisor, divisor);

		stack.push_back(_reciprocal);
		stack.push_back(_dividend);
		//dumpStack(1);
		Mul();
		//dumpStack(2);
		if (stack.back()->number.size())  for (int i = 0; i < shift;i++)
			Div2E30(*stack.back());

		while (1)
		{
			/*counter++;
			std::cout << "counter " << counter << std::endl;*/
			//if (stack.size() == 0)
			//	std::cout << "Something is rotten" << std::endl;
			/*if (stack.back()->number.size() > 3)
				std::cout << "Something is rotten" << std::endl;*/

			Dup(*Quotient, *stack.back());

			stack.push_back(_divisor);
			//dumpStack(3);
			Mul();
			//dumpStack(4);
			ChangeSign();
			stack.push_back(_dividend);
			Add();
			if (IsAbiggerNummerically(divisor, *stack.back()))
				break;

			stack.push_back(_reciprocal);
			//dumpStack(5);
			Mul();
			//dumpStack(6);
			for (int i = 0; i < shift;i++) {
				Div2E30(*stack.back());
			}
			if (IsZero(*stack.back())) {
				stack.pop_back();
				Push(1);
			};
			stack.push_back(Quotient);
			Add();
		}
		/* we are done */
		if (stack.back()->sign == -1)
		{
			stack.push_back(_divisor);
			Add();
			Push(-1); // adjust Quotient
			stack.push_back(Quotient);
			Add();
		}
		else
		{
			stack.push_back(Quotient);
		}
		Swap();
	}
}

void Calculator2E30::Div2E30(BInt2E30& A)
{
	switch (A.number.size())
	{
	case 0:
		break;

	case 1:
		A.number.pop_back();
		A.number.push_back(0);
		A.sign = 1;
		break;

	default:
		for (int ix = 0; A.number.size() && (ix < (A.number.size() - 1)); ix++)
		{
			A.number[ix] = A.number[ix + 1];
		}
		A.number.pop_back();
		break;
	}
}

void Calculator2E30::GCDAux(BInt2E30Ptr X, BInt2E30Ptr Y)
{
	BInt2E30 g;
	BInt2E30 x;  		Dup(x, *X);
	BInt2E30 y;  		Dup(y, *Y);

	/* count common 2^30 factors */
	int i = 0;
	u64 min_sz = std::min(x.number.size(), y.number.size());
	for (u64 i = 0; i < min_sz; i++)
		if ((x.number[i] | y.number[i]) == 0)  i++; /* counter #1000 */
		else break;

	/* then remove them and update g */
	for (int i2 = 0; i2 < i; i2++)  /* then remove them and update g */
	{
		g.number.push_back(0);
		Div2E30(x); Div2E30(y);
	}
	g.number.push_back(1);


	/* copy to temps*/
	BInt2E30 u; Dup(u, x);
	BInt2E30 v; Dup(v, y);

	/* remove common 2 factors */
	while ((u.number[0] & 1) == 0 && (v.number[0] & 1) == 0)
	{
		Div2(u); Div2(v); Mul2(g);
	};

	BInt2E30 A; A.number.push_back(1);
	BInt2E30 B; B.number.push_back(0);
	BInt2E30 C; C.number.push_back(0);
	BInt2E30 D; D.number.push_back(1);

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


	BInt2E30Ptr a1(new BInt2E30); Dup(*a1, C);
	BInt2E30Ptr b1(new BInt2E30); Dup(*b1, D);
	BInt2E30Ptr gcd(new BInt2E30); Dup(*gcd, v);
	BInt2E30Ptr gmul(new BInt2E30); Dup(*gmul, g);

	stack.push_back(a1);
	stack.push_back(b1);
	stack.push_back(gcd);
	stack.push_back(gmul);
	Mul();
}

void Calculator2E30::GCD() {
	if (stack.size() < 2)
		std::cout << "GCD() needs two arguments" << std::endl;
	else {
		if (IsZero(*stack[stack.size() - 1]))
			std::cout << "GCD() needs non-zero arguments" << std::endl;
		else if (IsZero(*stack[stack.size() - 2]))
			std::cout << "GCD() needs non-zero arguments" << std::endl;
		else {
			BInt2E30Ptr A = stack.back(); stack.pop_back();
			BInt2E30Ptr B = stack.back(); stack.pop_back();
			GCDAux(A, B);
		}
	}


}
void Calculator2E30::ChangeSign() {
	if (stack.size() > 0) {
		BInt2E30Ptr temp(new BInt2E30); Pop(*temp);
		temp->sign = -1 * temp->sign;
		stack.push_back(temp);
	}

}

/*
  The Asign and Bsign parameters are normally A and B real signs, but sometimes we
  want to overrule that.
*/

void Calculator2E30::AddAux( int Asign, BInt2E30& A, int Bsign, const BInt2E30& B)
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
			if (A.number[i] & RMODMASK) {
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

void Calculator2E30::Add() {
	if (stack.size() < 2)
		std::cout << "Add() needs two arguments" << std::endl;
	else {
		BInt2E30Ptr temp(new BInt2E30); 

		Dup( *temp, *stack.back());
		stack.pop_back();
		AddAux(temp->sign, *temp, stack.back()->sign, *stack.back());
		stack.pop_back();
		Normalize(temp);
		stack.push_back(temp);
	}
}

void Calculator2E30::Jacobi() {
	if (stack.size() < 2)
		std::cout << "Jacoby/Legendre  needs two arguments" << std::endl;
	else {
		BInt2E30Ptr Res(new BInt2E30); Res->number.clear();
		BInt2E30Ptr A(new BInt2E30); Pop(*A);
		BInt2E30Ptr M(new BInt2E30); Pop(*M);
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
					BInt2E30 temp;
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
void Calculator2E30::Div2(unsigned int Power) {

	unsigned int _Shift = Power;
	if (stack.size() > 0) {
		int Ctemp = 0;
		int C1 = 0;
		BInt2E30Ptr temp(new BInt2E30);
		Pop(*temp);

		while (_Shift > 30) //
		{
			for (int i = 1 ; i < (int)temp->number.size() ; i++)
				temp->number[i - 1] = temp->number[i];
			_Shift -= 30;
			if (temp->number.back() == 0) temp->number.pop_back();
		}
		while (_Shift > 0)
		{
			Ctemp = 0;
			C1 = 0;
			for (int i = (int)temp->number.size() - 1; i >= 0; i--)
			{
				Ctemp = temp->number[i] % 2;
				temp->number[i] = temp->number[i] / 2;
				temp->number[i] = (C1 * (RMOD / 2)) + temp->number[i];
				C1 = Ctemp;
			}
			_Shift--;
			if (temp->number.back() == 0) temp->number.pop_back();
		}
		Normalize(temp);

		stack.push_back(temp);
	}


}

#define BITS 30
void Calculator2E30::Rand() {
	if (stack.size() > 0) {

		int sz = BITS * (uint)(stack.back()->number.size() - 1);
		int s = stack.back()->number.back();

		if (s == 0) {
			std::cout << "Zero Argument " << std::endl;
			return;
		}

		while (s > 0) { s = s>>1 ; sz++; } //count the numbers of digits in the MSInt

		while (sz > 0) {
			// approx 90% of all numbers smaller than sz is  in range sz ... sz/10
			// approx 99% of all numbers smaller than sz is  in range sz ... sz/100
			// approx 99.9% of all numbers smaller than sz is  in range sz ... sz/1000.....
			if (_Rand(RMOD) >= (RMOD /2)) break;
			sz--;
		}
		int size = std::max(1, sz);

		BInt2E30Ptr temp(new BInt2E30); temp->number.clear();
		for (; size > BITS; size = size - BITS)
			temp->number.push_back(_Rand(RMOD)); //  _Rand() returns an integer in the range 0..RMOD-1


		int i = _Rand(RMOD);

		while (BITS > size) {
			i = i / 2;
			size++;
		}
		i = i % stack.back()->number.back();

		temp->number.push_back(i);

		Normalize(temp);
		if (stack.back()->number.size() <= temp->number.size())
		{
			if (stack.back()->number.size() < temp->number.size())
				std::cout << "Rand() something is rotten " << std::endl;
			// debug trap
			else if (stack.back()->number[stack.back()->number.size() - 1] <= temp->number[temp->number.size() - 1])
				std::cout << "Rand() something is rotten " << std::endl;
		}
		stack.pop_back();
		stack.push_back(temp);
	}
	else
		std::cout << "Empty Stack !" << std::endl;
}


// 
/* store operations  */
void Calculator2E30::PopStore(const std::string& loc) {
	if (stack.size() > 0) {
		u64 t = stack.back()->number.size();  // debug

		Store[loc] = stack.back();
		stack.pop_back();
		if (t != Store[loc]->number.size())  // debug
			std::cout << "PopStore Error" << std::endl;  // debug
	}

}
void Calculator2E30::PushStore(const std::string& loc) {
	if (Store.find(loc) != Store.end())
		stack.push_back(Store.find(loc)->second);
	else
		std::cout << "Lookup failed" << std::endl;
}
void Calculator2E30::ClearStore(const std::string& loc) {
	if (Store.find(loc) != Store.end())
		Store.erase(loc);
	else
		std::cout << "Lookup failed" << std::endl;
}

/*  Various Predicates */
int  Calculator2E30::TOSStatus() {
	if (stack.size() > 0)
		if (IsZero(*(stack.back()))) return 0;
		else return stack.back()->sign;
	return 0;

}
int  Calculator2E30::TOSSize() {
	if (stack.size() > 0) return (int)stack.back()->number.size();
	return 0;
}
bool Calculator2E30::IsLarger() {
	size_t sz = stack.size();
	if (sz < 2)
		std::cout << "IsLarger(): not enough arguments" << std::endl;
	else
	{
		if ((stack[sz - 1]->sign == 1) && (stack[sz - 2]->sign == -1)) return true;
		if ((stack[sz - 1]->sign == 1) && (stack[sz - 2]->sign == 1))
			return IsAbiggerNummerically(stack[sz - 1], stack[sz - 2]);
		if ((stack[sz - 1]->sign == -1) && (stack[sz - 2]->sign == -1))
			return IsAbiggerNummerically(stack[sz - 2], stack[sz - 1]);
	}
	return false;
}
bool Calculator2E30::IsEqual() {
	size_t sz = stack.size();
	if (sz < 2)
	{
		std::cout << "IsLarger(): not enough arguments" << std::endl;
		return false;
	}
	else
		return IsEqual(*stack[sz - 1], *stack[sz - 2]);
	return true;
}

bool Calculator2E30::IsEqual(BInt2E30 A, BInt2E30 B)
{
	if (A.sign != B.sign) return false;
	if (A.number.size() != B.number.size()) return false;
	for (size_t i = 0; i < A.number.size(); i++)
		if (A.number[i] != B.number[i]) return false;
	return true;
}


bool Calculator2E30::IsZero() {
	if (stack.size() > 0) 	 return stack.back()->sign == 1
		&& stack.back()->number[0] == 0
		&& stack.back()->number.size() == 1;
	return false;
}

bool Calculator2E30::IsOne() {
	if (stack.size() > 0) return stack.back()->sign == 1
		&& stack.back()->number[0] == 1
		&& stack.back()->number.size() == 1;
	return false;
}
bool Calculator2E30::IsMinusOne() {
	if (stack.size() > 0) 	return stack.back()->sign == -1
		&& stack.back()->number[0] == 1
		&& stack.back()->number.size() == 1;
	return false;

}
bool Calculator2E30::IsEven() {
	if (stack.size() > 0) 	return ((stack.back()->number[0] & 1) == 0);
	return false;
}


// for internal use......
void Calculator2E30::dumpStack(int p) {
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


bool Calculator2E30::IsZero(BInt2E30 arg) {
	int num = 0;
	for (unsigned int i = 0; !num && (i < arg.number.size()); i++)	num |= arg.number[i];
	return num == 0;

}


bool Calculator2E30::IsAbiggerNummerically(BInt2E30Ptr A, BInt2E30Ptr B) {
	return IsAbiggerNummerically(*A, *B);
}

bool Calculator2E30::IsAbiggerNummerically(BInt2E30 A, BInt2E30 B)
{
	if (IsZero(A) && IsZero(B)) return false;
	if (IsZero(A) && !IsZero(B)) return false;
	if (!IsZero(A) && IsZero(B)) return true;
	/* they are both positive */
	if (A.number.size() > B.number.size()) return true;
	if (A.number.size() < B.number.size()) return false;
	/* they are equal in size */
	for (u64 i = A.number.size(); i > 0; i--) {
		if (A.number[i - 1] > B.number[i - 1]) return true;
		if (A.number[i - 1] < B.number[i - 1]) return false;
	}
	return false;
}

bool Calculator2E30::IsALarger(BInt2E30 A, BInt2E30 B) {

	if (IsAbiggerNummerically(A, B))
		if (A.sign == 1) return true;
		else return false;
	else
		if (A.sign == -1) return true;
		else return false;
	return false;
}



void Calculator2E30::Div2(BInt2E30& A) {

	int borrow = 0;
	for (s64 ix = A.number.size() - 1; ix >= 0; ix--)
	{
		int t = A.number[ix] + borrow; borrow = 0;
		if (t & 1)  borrow = RMOD;
		A.number[ix] = t >> 1;
	}
	while (A.number.size() && A.number.back() == 0) A.number.pop_back();
}

void Calculator2E30::Mul2(BInt2E30& A) {
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

	//for (u64 i = 0; i < A.number.size(); i++) A.number[i] = A.number[i] + A.number[i];

	//for (u64 i = 0; i < A.number.size(); i++)
	//{
	//	A.number[i] += carry;
	//	if (A.number[i] & RMODMASK) {
	//		A.number[i] -= RMOD;
	//		carry = 1;
	//	}
	//	else
	//		carry = 0;
	//}
	if (carry) A.number.push_back(carry);
}
void Calculator2E30::Dup(BInt2E30& D, const BInt2E30& S) {

	D.number.clear();
	D.sign = S.sign;
	for (int ix = 0; ix < S.number.size(); ix++)
		D.number.push_back(S.number[ix]);
}

void Calculator2E30::Normalize(BInt2E30& b) {
	while (b.number.size() && b.number.back() == 0) b.number.pop_back();
	if (b.number.size() == 0) {
		b.number.push_back(0);
		b.sign = 1;
	}
}

#define FORMATSTRING "%08X"
void Calculator2E30::DumpInt(std::string name, const BInt2E30& arg) {
	char buffer[12];
	std::string* s = new std::string();
	for (int i = 0; i < arg.number.size(); i++) {
		sprintf(buffer, FORMATSTRING, arg.number[i]);
		char* c1 = buffer;
		char* c2 = buffer + 7 ;
		while (c1 < c2)
		{
			char t = *c1;
			*c1 = *c2;
			*c2 = t;
			c1++; c2--;
		}
		s->append(buffer);
	}
	while (s->size() && (s->back() == '0')) s->pop_back();
	if (s->size() == 0) s->append("0");
	if (arg.sign == -1) s->append("-");
	std::reverse(s->begin(), s->end());
	std::cout << name << " : " << *s << " ( " << s->length() << " ) " << std::endl;

}

#define SIMPLEMULT(s,A,sb, B)  SimpleAdditionSubtractionLadder1(s,A,sb, B)

void Calculator2E30::RussianPeasantMult() {
	BInt2E30 x;
	BInt2E30 y;

	Pop(x); 	Pop(y);
#if TESTMUL == 1
	std::cout << "RussianPeasantMult() " << std::endl;
	DumpInt(" x: ", x);
	DumpInt(" y: ", y);
#endif
	if (IsZero(x) || IsZero(y))
	{
		BInt2E30Ptr result(new BInt2E30());
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

void Calculator2E30::Exp()
{
	if (stack.size() < 2)
		std::cout << "Exp() needs two arguments" << std::endl;
	else {
		BInt2E30     Exponent;  Pop(Exponent);
		BInt2E30Ptr  Argument(new BInt2E30); Pop(*Argument);
		BInt2E30Ptr  Result(new BInt2E30); Result->number.push_back(1);
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
void Calculator2E30::SimpleAdditionSubtractionLadder1(int sign, s64 A, int BSign, const BInt2E30& B) {
{
		/* this is taken from Crandall& Pomerance
		*   "Prime Numbers,  A Computational Perspective" 2nd edition
		*/
		BInt2E30Ptr result(new BInt2E30); result->number.clear();result->number.push_back(0);
		

#if TESTMUL == 1
		std::cout << "SimpleAdditionSubraction: " << std::endl;
		std::cout << "Argument A:  " << A << std::endl;
		DumpInt("Argument B:  ", B);
#endif

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
					AddAux(1, *result,1, B);
				}
				else {
					AddAux(1, *result, -1, B);
				}
			}
			Mask = Mask >> 1;
		}
		Normalize(*result);
		if (!IsZero(*result)) result->sign = sign;

#if TESTMUL == 1
		DumpInt("result:  ", *result);
#endif
		stack.push_back(result);
	}

}


uint  Calculator2E30::_Rand(uint UpperBound) {
	uint ix2 = dist->operator()(rd);
	return ix2 % UpperBound;
}

std::string* Calculator2E30::ItoA()
{
	Calculator c;
	BInt  t10;
	BInt2E30 t2;
	Pop(t2);
	Convert2E30to10E9(t10, t2);
	c.Push(t10);
	return c.ItoA();

}

