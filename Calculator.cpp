/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/

#include <ctype.h>
#include "Calculator.h"

/* radix definitions  */
#define RMOD 1000000000
#define RMOD3 1000
#define FORMATSTRING "%09d"
#define DIGITS 9

Calculator::Calculator(){
	dist = new std::uniform_int_distribution<uint>(0, RMOD - 1);
}

Calculator::~Calculator() {
	delete dist;
}

int Calculator::TOSStatus() {
	if (stack.size() > 0)
		if(IsZero(*(stack.back()))) return 0;
		else return stack.back()->sign;
	return 0;
}

int Calculator::TOSSize() {
	if (stack.size() > 0) return (int) stack.back()->number.size();
	return 0;
}

bool Calculator::IsZero(BInt arg)
{
	int num = 0;
	for (unsigned int i = 0; !num && ( i < arg.number.size()); i++)	num |= arg.number[i];
	return num == 0;
}

bool Calculator::IsZero()
{
	if (stack.size() > 0) 	 return stack.back()->sign == 1
		&& stack.back()->number[0] == 0
		&& stack.back()->number.size() == 1;
	return false;
}

bool Calculator::IsOne()
{
	if (stack.size() > 0) return stack.back()->sign == 1 
		 && stack.back()->number[0]==1 
		 && stack.back()->number.size() == 1;
	return false;
}

bool Calculator::IsMinusOne()
{
	if (stack.size() > 0) 	return stack.back()->sign == -1
		&& stack.back()->number[0] == 1
		&& stack.back()->number.size() == 1;
	return false;
}


bool Calculator::IsEven()
{
	if (stack.size() > 0) 	return ((stack.back()->number[0] & 1) == 0);
	return false;
}

bool Calculator::IsLarger() {
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

bool Calculator::IsEqual() {
	size_t sz = stack.size();
	if (sz < 2)
	{
		std::cout << "IsLarger(): not enough arguments" << std::endl;
		return false;
	}
	else
	{
		if ((stack[sz - 1]->sign != stack[sz - 2]->sign)) return false;
		if ((stack[sz - 1]->number.size() != stack[sz - 2]->number.size())) return false;
		for(size_t i = 0 ; i < stack[sz - 1]->number.size(); i++)
			if ((stack[sz - 1]->number[i] != stack[sz - 2]->number[i])) return false;
	}
	return true;
}

bool Calculator::IsEqual(BInt A, BInt B) 
{
	if (A.sign != B.sign) return false;
	if (A.number.size() != B.number.size()) return false;
	for (size_t i = 0; i < A.number.size(); i++)
		if (A.number[i] != B.number[i]) return false;
	return true;
}

bool Calculator::IsALarger(BInt A, BInt B)
{
	if( IsAbiggerNummerically(A,B))
		if (A.sign == 1) return true;
		else return false;
	else
		if (A.sign == -1) return true;
		else return false;
	return false;
}

bool Calculator::IsAbiggerNummerically(BIntPtr A, BIntPtr B) { return IsAbiggerNummerically(*A, *B); }

bool Calculator::IsAbiggerNummerically(BInt A, BInt B)
{
	if (IsZero(A) && IsZero(B)) return false;
	if (IsZero(A) && !IsZero(B)) return false;
	if (!IsZero(A) && IsZero(B)) return true;
	/* they are both positive */
	if (A.number.size() > B.number.size()) return true;
	if (A.number.size() < B.number.size()) return false;
	/* they are equal in size */
	for (u64 i = A.number.size(); i > 0; i--) {
		if (A.number[i-1] > B.number[i-1]) return true;
		if (A.number[i-1] < B.number[i-1]) return false;
	}
	return false;
}

uint  Calculator::_Rand(uint UpperBound) 
{
	uint ix2 = dist->operator()(rd);
	return ix2 % UpperBound;
}

void Calculator::Rand(){  // we are not aiming for perfection here.....
	if (stack.size() > 0) {

		int sz =  DIGITS * (uint)( stack.back()->number.size()-1);
		int s = stack.back()->number.back();

		if (s == 0) {
			std::cout << "Zero Argument " << std::endl;
			return;
		}

		while (s > 0) { s = s / 10; sz++; } //count the numbers of digits in the MSInt

		while (sz > 0) {  
			// approx 90% of all numbers smaller than sz is  in range sz ... sz/10
			// approx 99% of all numbers smaller than sz is  in range sz ... sz/100
			// approx 99.9% of all numbers smaller than sz is  in range sz ... sz/1000.....
			if (_Rand(RMOD) >= (RMOD / 10)) break;
			sz--;
		}
		int size = std::max(1,sz); 

		BIntPtr temp(new BInt); temp->number.clear();
		for (; size >= DIGITS; size = size - DIGITS)
			temp->number.push_back(_Rand(RMOD)); //  _Rand() returns an integer in the range 0..RMOD-1


		int i =  _Rand(RMOD); 

		while (DIGITS > size ) {
			i = i / 10;
			size++;
		}
		i = i% stack.back()->number.back();

		temp->number.push_back(i);

		Normalize(temp);
		stack.pop_back();
		stack.push_back(temp);
	}
	else 
		std::cout << "Empty Stack !" << std::endl;

}

/*inline*/ void Calculator::PopStore(const std::string& loc) 
{ 

	if (stack.size() > 0){ 
		u64 t = stack.back()->number.size();  // debug

		Store[loc] = stack.back();	
		stack.pop_back(); 
		if (t != Store[loc]->number.size())  // debug
			std::cout << "PopStore Error" << std::endl;  // debug
	}
};


void Calculator::Div2(unsigned int Power)
{
	unsigned int _Shift = Power;
	if (stack.size() > 0) {
		int Ctemp = 0;
		int C1 = 0;
		BIntPtr temp(new BInt);
		Pop(*temp);

		while (_Shift > 9) // RMOD == 2^9 * 5^9
		{
			Ctemp = 0;
			C1 = 0;
			for (int i = (int) temp->number.size() - 1; i >= 0; i--)
			{
				Ctemp = temp->number[i] % 512;
				temp->number[i] = temp->number[i] / 512;
				temp->number[i] = (C1 * (RMOD / 512)) + temp->number[i];
				C1 = Ctemp;
			}
			_Shift -= 9;
			if (temp->number.back() == 0) temp->number.pop_back();
		}
		while (_Shift > 0)
		{
			Ctemp = 0;
			C1 = 0;
			for (int i = (int) temp->number.size() - 1; i >= 0; i--)
			{
				Ctemp = temp->number[i] % 2;
				temp->number[i] = temp->number[i] / 2;
				temp->number[i] = (C1 * (RMOD /2)) + temp->number[i];
				C1 = Ctemp;
			}
			_Shift--;
			if (temp->number.back() == 0) temp->number.pop_back();
		}
		Normalize(temp);

		stack.push_back(temp);
	}
}


void Calculator::Push(const BInt  &b)
{
	BIntPtr temp(new BInt);
	Dup(*temp, b);
	stack.push_back(temp);
}


void Calculator::Push(int  i)
{
	BIntPtr temp(new BInt); temp->number.clear();
	s64 it = i;
	temp->sign = (it >= 0)? 1 : - 1;

	if (it < 0) it = -1 * i;
	do {
		temp->number.push_back(it % (s64)RMOD);
		it = it / RMOD;
	} while (it > 0);
	stack.push_back(temp);
}

void Calculator::Push( char *c)
{	
	BIntPtr temp(new BInt);
	char *ct1 , *ct2 ;

	ct1 = ct2 = c ;
	temp->sign = 1;

	while (isspace(*ct1)) ct1++;
	if (*ct1 == '-') {
		temp->sign = -1;
		ct1++;
	}
	ct2 = ct1;

	while (isdigit(*ct2)) ct2++; // find the end;
	if (ct2 == ct1) {  //bail out something is wrong
		std::cerr << "unknown format " << c << std::endl;
	}
	// we assume we have something number like 
	temp->number.clear();
	int tempInt = 0;
	while ((ct2 -ct1) > DIGITS -1) {
		tempInt = 0;
		ct2 = ct2 - DIGITS;
		for (int i = 0; i < DIGITS; i++) tempInt = tempInt * 10 + ct2[i] - '0';
		temp->number.push_back(tempInt);
	}
	tempInt = 0;
	for( int i = 0; i < ct2-ct1; i++) tempInt = tempInt * 10 + ct1[i] - '0';
	temp->number.push_back(tempInt);

	Normalize(temp);
	stack.push_back(temp);
}

std::string* Calculator::ItoA()
{
	char buffer[12];
	std::string* s = new std::string();
	if (stack.size() > 0) {
		BIntPtr o = stack.back(); stack.pop_back();
		for (int i = 0; i < o->number.size(); i++) {
			sprintf(buffer, FORMATSTRING, o->number[i]);
			char* c1 = buffer;
			char* c2 = buffer + DIGITS-1;
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
		if (o->sign == -1) s->append("-");
		std::reverse(s->begin(), s->end());
		return s;
	}
	s->append("Empty Stack");
	return s;
}

void Calculator::ChangeSign() 
{
	if (stack.size() > 0) {
		BIntPtr temp(new BInt); Pop(*temp);
		temp->sign = -1 * temp->sign;
		stack.push_back(temp);
	}
}

/*

the code below is adapted from

		Algorithm 14.9 	in

		The Handbook of Applied Cryptograpy by A.Menezes, P van Oorschot and S.Vanstone (CRC Press 1996).
*/
void Calculator::Add()
	{
		if (stack.size() < 2)
			std::cout << "Add() needs two arguments" << std::endl;
		else {
			/* setup  */
			u64 min_sz = std::min(stack[stack.size() - 1]->number.size(), stack[stack.size() - 2]->number.size());
			BIntPtr temp(new BInt); temp->number.clear();

			if (stack[stack.size() - 1]->sign == stack[stack.size() - 2]->sign)
			{
				/* they have identical signs */
				/* this is done in three steps
				*    1)  copy TOS to temp
				*    2)  Add the digits without carry
				*    3)  do ' + carry mod RMOD'
				*/

				temp->sign = stack[stack.size() - 1]->sign;
				/* step 1 */
				for (uint i = 0; i < stack.back()->number.size(); i++) temp->number.push_back(  stack.back()->number[i]);

				stack.pop_back();

				/* step 2: add the other number, it may be longer or shorter than temp */
				for (u64 i = 0; i < min_sz; i++)  temp->number[i] = temp->number[i] + stack.back()->number[i];

				if (stack.back()->number.size() > temp->number.size()) 
					for (u64 i = min_sz; i < stack.back()->number.size(); i++)
						temp->number.push_back(stack.back()->number[i]);

				stack.pop_back(); // done with all inputs

				/* step 3 */
				int carry = 0;
				for (u64 i = 0; i < temp->number.size(); i++)
				{
					temp->number[i] += carry;
					if (temp->number[i] >= RMOD) {
						temp->number[i] -= RMOD;
						carry = 1;
					}
					else
						carry = 0;
				}
				if (carry) temp->number.push_back(carry);

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
				if (stack.back()->sign == 1) Swap();  // we put the negative number on top

				/* step 1: negate digits in the negative number, which is TOS,  and save them in temp  */
				for (uint i = 0; i < stack.back()->number.size(); i++) temp->number.push_back( - stack.back()->number[i]);

				stack.pop_back(); // remove the negative number from top

				/* step 2: add the other number, which is positive, it may be longer or shorter than temp */
				for (u64 i = 0; i < min_sz; i++)  temp->number[i] = temp->number[i] + stack.back()->number[i];

				if (stack.back()->number.size() > temp->number.size())  /* is the positive number bigger than temp */
					for (u64 i = min_sz; i < stack.back()->number.size(); i++)	temp->number.push_back(stack.back()->number[i]);

				stack.pop_back(); // done with all inputs

				/* step 3: now adjust carry/borrow and mod */
				int c = 0;
				for (int i = 0; i < temp->number.size(); i++)
				{
					temp->number[i] += c;
					if (temp->number[i] < 0) {
						c = -1;
						temp->number[i] += RMOD;
					}
					else
						c = 0;
				}

				/* step 4: final correction if c is set */
				if (c != 0) {
					temp->sign = -1;
					int c1 = 0;
					for (int i = 0; i < temp->number.size(); i++)
					{
						temp->number[i] = 0 - temp->number[i] + c1;
						if (temp->number[i] < 0) {
							c1 = -1;
							temp->number[i] += RMOD;
						}
						else
							c1 = 0;
					}
				}
			}
			Normalize(temp);
			stack.push_back(temp);
		}
}

void Calculator::PushStore(const std::string& loc) 
{
	if( Store.find(loc) != Store.end())
			stack.push_back( Store.find(loc)->second);
	else
		std::cout << "Lookup failed" << std::endl;
}

void Calculator::ClearStore(const std::string& loc)
{
	if (Store.find(loc) != Store.end())
		Store.erase(loc);
	else
		std::cout << "Lookup failed" << std::endl;
}

void Calculator::Swap()
{

	if (stack.size() > 1) {

		BIntPtr t1 = stack.back(); 		stack.pop_back();
		BIntPtr t2 = stack.back();		stack.pop_back();
		stack.push_back(t1);
		stack.push_back(t2);
	}
}

#define OVERALLOCATION 2

void Calculator::Mul()
{

	if (stack.size() < 2)
		std::cout << "Mul() needs two arguments" << std::endl;
	else if (SMALLNUMBERLIMIT && stack[stack.size() - 1]->number.size() < 1+SMALLNUMBERLIMIT)
		RussianPeasantMult();
	else if (SMALLNUMBERLIMIT && stack[stack.size() - 2]->number.size() < 1+SMALLNUMBERLIMIT)
		RussianPeasantMult();
	else {
		PrimeFactorDFT pf;
		BIntPtr temp(new BInt); temp->number.clear();

		temp->sign = stack[stack.size() - 1]->sign * stack[stack.size() - 2]->sign;

		factorSeq  factors;

		u64 min_sz = stack[stack.size() - 1]->number.size() +
			stack[stack.size() - 2]->number.size();

		u64 Length = 1;
		
		int factorlist[] = { 31,19,17,13,11,7,5,3,2 };

		for( int i = 0; i < sizeof(factorlist)/sizeof(factorlist[0]); i++)
		{ 
			factors.push_back(factorlist[i]); 
			Length *= factors.back();
			if (Length > 3* 2 * min_sz) break; /* 3 because we go from Radix 10^9 to Radix 10^3 */
			if (i == (sizeof(factorlist) / sizeof(factorlist[0]) - 1)) {
				std::cout << "Sorry Number too Large" << std::endl;
				return;
			}
		}

		pf.SetFactors(factors);

		if (pf.Status() > 0) {
			/* this should be less dynamic, but right now it is OK */
			Data* real1 = new Data[pf.Status()+OVERALLOCATION];
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

			real3[0] = real1[0]* imag1[0];
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

			/* convert back to radix 10^9 from radix 10^3 double*/

			for (int i = 0; i < pf.Status(); i+=3)
			{
				int t0 = (int)real1[i];
				int t1 = 1000 * (int)real1[i+1];  // these values are 0 when we reach
				int t2 = 1000*1000* (int)real1[i+2];// the end of buffer, due to
				temp->number.push_back(t0+t1+t2);//OVERALLOCATION
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

void Calculator::LoadFFT(BIntPtr A, Data* Buffer)
{
	int carry = 0;
	int FFTIndex = 0;
	for (int ix = 0; ix < A->number.size(); ix++) {
		int temp = A->number[ix];   // temp is  radix 10^9 
		for (int i = 0; i < 3; i++) {
			int tmp = (temp%RMOD3);// tmp is  radix 10^3
			temp = temp / RMOD3;
			tmp = tmp + carry; carry = 0;
			if (tmp > (RMOD3 / 2) - 1) {  // tmp is 'balanced' radix 10^3 
				tmp = tmp - RMOD3;
				carry++;
			}
			Buffer[FFTIndex++] = (double)tmp;
		}
	}
	if (carry)
		Buffer[FFTIndex] = 1.0;
}

void Calculator::Carry(s64 size, Data* Buffer)
{
	s64 carry = 0;
	s64 tmp = 0;

	/* conversion from 'balanced' Radix 10^3 double  to unbalanced Radix 10^3 double */
	for (s64 i = 0; i < size; i++)
	{
		tmp = (s64) std::round(Buffer[i]);

		tmp += carry; carry = 0;
		while (tmp < ( - 1 * (RMOD3 / 2))) { tmp += RMOD3; carry--;	}
		while (tmp > ((RMOD3 / 2)-1)) { 	tmp -= RMOD3; carry++;	}
		Buffer[i] = (double) tmp;
	}
	carry = 0;
	/* carry we are still in Radix 10^3 double*/
	for (s64 i = 0; i < size; i++)
	{
	    tmp = (s64) std::round(Buffer[i]);
		tmp += carry;  carry = 0;
		if (tmp < 0) { tmp += RMOD3; carry--; }
		Buffer[i] = double(tmp) ;
	}

}

void Calculator::GCD()
{
	if (stack.size() < 2)
		std::cout << "GCD() needs two arguments" << std::endl;
	else {
		if(IsZero( *stack[stack.size()-1])) 
			std::cout << "GCD() needs non-zero arguments" << std::endl;
		else if (IsZero(*stack[stack.size() - 2]))
			std::cout << "GCD() needs non-zero arguments" << std::endl;
		else {
			BIntPtr A = stack.back(); stack.pop_back();
			BIntPtr B = stack.back(); stack.pop_back();
			GCDAux(A, B);
		}
	}
}

void Calculator::Div2(BInt &A)
{
	int borrow = 0;
	for (s64 ix = A.number.size() - 1; ix >= 0; ix--)
	{
		int t = A.number[ix] + borrow; borrow = 0;
		if (t & 1)  borrow = RMOD;
		A.number[ix] = t >> 1;
	}
	while (A.number.size() && A.number.back() == 0) A.number.pop_back();

}

void Calculator::Div1000000000(BInt& A)
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


void Calculator::Mul2(BInt &A)
{
	int carry = 0; 
	for (int ix = 0; ix < A.number.size(); ix++)
	{
		int t = (A.number[ix]  << 1) + carry;
		carry = 0;
		if (t >= RMOD) {
			carry = 1;
			t -= RMOD;
		}
		A.number[ix] = t;
	}
	if (carry) A.number.push_back(1);
}

void Calculator::Dup(BInt& D, const BInt& S)
{
	D.number.clear();
	D.sign = S.sign;
	for (int ix = 0; ix < S.number.size(); ix++)
		D.number.push_back(S.number[ix]);
}

void Calculator::GCDAux(BIntPtr X, BIntPtr Y)
{
	BInt g;
	BInt x;  		Dup(x, *X);
	BInt y;  		Dup(y, *Y);

	/* remove common 1000 factors */
	int i = 0;
	if (x.number.size() > 1 && y.number.size() > 1)
		while ((x.number[i] | y.number[i]) == 0)  i++; /* counter #1000 */

	/* then remove them and update g */
	for (int i2 = 0; i2 < i; i2++)  /* then remove them and update g */
	{
		g.number.push_back(0);
		Div1000000000(x); Div1000000000(y);
	}
	g.number.push_back(1);


	/* copy to temps*/
	BInt u; Dup(u, x);
	BInt v; Dup(v, y);

	/* remove common 2 factors */
	while((u.number[0] & 1) == 0 && (v.number[0] & 1) == 0)
		{
			Div2(u); Div2(v); Mul2(g);
		};

	BInt A; A.number.push_back(1);
	BInt B; B.number.push_back(0);
	BInt C; C.number.push_back(0);
	BInt D; D.number.push_back(1);

	while (!IsZero(u)){
		//std::cout << "-----------------------------------" << std::endl;
		DUMPINT("u", u); 		DUMPINT("v", v);		DUMPINT("A", A);		
		DUMPINT("B", B);		DUMPINT("C", C);		DUMPINT("D", D);
		//std::cout << "-----------------------------------" << std::endl;

		while ((u.number[0] & 1) == 0)
		{
			Div2(u); DUMPINT("u", u);
			if ((A.number[0] & 1) == 0 && (B.number[0] & 1) == 0) {	DUMPINT("A", A);
				Div2(A); DUMPINT("B", B);
				Div2(B);
			}
			else {				DUMPINT("A", A);
				Add(A, y);		DUMPINT("A", A);
				Div2(A);		DUMPINT("A", A);DUMPINT("B", B);
				Sub(B, x);		DUMPINT("B", B);
				Div2(B);		DUMPINT("B", B);
			}
		}

		while ((v.number[0] & 1) == 0)	{  DUMPINT("v", v);
			Div2(v); DUMPINT("v", v);
			if ((C.number[0] & 1) == 0 && (D.number[0] & 1) == 0) {	DUMPINT("C", C);
				Div2(C);  DUMPINT("D", D);
				Div2(D);
			}
			else {	DUMPINT("C", C);  DUMPINT("y", y);
				Add(C, y);	DUMPINT("C", C);
				Div2(C);	DUMPINT("C", C);  DUMPINT("D", D);	DUMPINT("x", x);
				Sub(D, x);	DUMPINT("D", D);
				Div2(D);	DUMPINT("D", D);
			}
		}
		DUMPINT("u", u);	DUMPINT("v", v);

		if (IsALarger(u, v) || IsEqual(u, v))
		{	DUMPINT("u", u);DUMPINT("v", v);
			Sub(u, v);	DUMPINT("u", u);	DUMPINT("A", A);	DUMPINT("C", C);
			Sub(A, C);	DUMPINT("A", A);	DUMPINT("B", B);	DUMPINT("D", D);
			Sub(B, D);	DUMPINT("B", B);
		}
		else {	DUMPINT("u", u);   DUMPINT("v", v);
			Sub(v, u);	DUMPINT("v", v);  DUMPINT("A", A);	DUMPINT("C", C);
			Sub(C, A);	DUMPINT("C", C);  DUMPINT("B", B);	DUMPINT("D", D);
			Sub(D, B);	DUMPINT("D", D);

		}
	};


	BIntPtr a1(new BInt); Dup(*a1, C);
	BIntPtr b1(new BInt); Dup(*b1, D);
	BIntPtr gcd(new BInt); Dup(*gcd, v);
	BIntPtr gmul(new BInt); Dup(*gmul, g);
	
	stack.push_back(a1);
	stack.push_back(b1);
	stack.push_back(gcd);
	stack.push_back(gmul);
	Mul();
}

void Calculator::DumpInt(std::string name,  const BInt& arg)
{
	char buffer[12];
	std::string* s = new std::string();
	for (int i = 0; i < arg.number.size(); i++) {
		sprintf(buffer, FORMATSTRING, arg.number[i]);
		char* c1 = buffer;
		char* c2 = buffer + DIGITS - 1;
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
	std::cout << name <<" : " << *s << " ( " << s->length() << " ) " << std::endl;

}

void Calculator::QuotientRemainder()
{
	BIntPtr Quotient(new BInt); Quotient->number.clear();
	Quotient->number.push_back(0);
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
		BIntPtr    Remainder(new BInt); Remainder->number.clear();
		Quotient->number.push_back(stack[stack.size() - 2]->number[0] / stack[stack.size() - 1]->number[0]);
		Remainder->number.push_back(stack[stack.size() - 2]->number[0] % stack[stack.size() - 1]->number[0]);
		stack.pop_back(); stack.pop_back();
		stack.push_back(Quotient);
		stack.push_back(Remainder);
	}
	else {

		//dumpStack(0);
		BInt    divisor;  Pop(divisor);
		BInt    dividend; Pop(dividend);
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
		int     reciprocal = RMOD / divisor.number.back();
		BIntPtr _reciprocal(new BInt); _reciprocal->number.clear(); _reciprocal->number.push_back(reciprocal);
		int     shift = (int)divisor.number.size();
		BIntPtr _dividend(new BInt); Dup(*_dividend, dividend);
		BIntPtr _divisor(new BInt);  Dup(*_divisor, divisor);

		stack.push_back(_reciprocal);
		stack.push_back(_dividend);
		//dumpStack(1);
		Mul(); 	
		//dumpStack(2);
		if (stack.back()->number.size())  for (int i = 0; i < shift;i++)
			Div1000000000(*stack.back());	

		while(1)
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
				Div1000000000(*stack.back()); 
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

#define SIMPLEMULT(S,A,B) RussianPeasantMultAux(S,A,B)

void Calculator::RussianPeasantMultAux(int sign, s64 A, const BInt& B)
{
	BIntPtr result(new BInt); result->number.clear();result->number.push_back(0);
	BInt    addend;  Dup(addend, B);
	if(!IsZero(B))
		while (A) {
			if (A & 1) Add(*result, addend);
			A = A >> 1;
			if (A) Mul2(addend);
		}
	else 
		result->number.push_back(0);
	Normalize(*result);
	if (!IsZero(*result)) result->sign = sign;
	else result->sign = 1; // 0 is positive.
	stack.push_back(result);
}

void Calculator::RussianPeasantMult()
{
	BInt x; 
	BInt y;
	
	Pop(x); 	Pop(y);
	s64 xint = 0;
	s64 yint = 0;

	if (x.number.size() < 1+ SMALLNUMBERLIMIT)
		for (s64 ix = x.number.size(); ix > 0;ix--)
			xint = (xint * RMOD) + x.number[ix - 1];
	if (y.number.size() < 1+ SMALLNUMBERLIMIT)
		for (s64 iy = y.number.size(); iy > 0;iy--)
			yint = (yint * RMOD) + y.number[iy - 1];

	if (yint == 0 )  		SIMPLEMULT(x.sign * y.sign, xint, y);
	else if (xint == 0)             SIMPLEMULT(x.sign * y.sign, yint, x);
	else if (xint < yint) 	        SIMPLEMULT(x.sign*y.sign, xint, y);
		  else            	SIMPLEMULT(x.sign * y.sign, yint, x);
}

void Calculator::Exp()
{
	if (stack.size() < 2)
		std::cout << "Exp() needs two arguments" << std::endl;
	else {
		BInt     Exponent;  Pop(Exponent);
		BIntPtr  Argument(new BInt); Pop(*Argument);
		BIntPtr  Result(new BInt); Result->number.push_back(1);
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

void Calculator::Jacobi()
{
	if (stack.size() < 2)
		std::cout << "Jacoby/Legendre  needs two arguments" << std::endl;
	else {
		BIntPtr A(new BInt); Pop(*A);
		BIntPtr M(new BInt); Pop(*M);
		BIntPtr Res(new BInt); Res->number.clear();

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
				BInt temp;
				Dup(temp, *A);
				Dup(*A, *M);
				Dup(*M, temp);
				if ((3 == (A->number[0] & 0x3)) && (3 ==( M->number[0] & 0x3)))
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

void Calculator::Normalize(BInt& b)
{
	while (b.number.size() && b.number.back() == 0) b.number.pop_back();
	if (b.number.size() == 0) {
		b.number.push_back(0);
		b.sign = 1;
	}
}



void Calculator::dumpStack(int p)
{
	printf("stack(%d): \n", p);
	if (stack.size() == 0) 
		std::cout << "empty stack" << std::endl;
	else
		for (int ix = 1; ix <= stack.size();ix++) {
			printf("level %-2d : ", ix-1);
			DUMPINT("", *stack[stack.size() - ix]);
		}
	printf("--------------------------\n");
}
