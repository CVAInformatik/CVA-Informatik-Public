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
//
#define CALCULATOR Calculator 
#define BINT BInt
#define BINTPTR BIntPtr
#define DIVRMOD(x) Div1000000000(x)

Calculator::Calculator() {
	stack.reserve(10);
	dist = new std::uniform_int_distribution<uint>(0, RMOD - 1);
}

Calculator::~Calculator() {
	delete dist;
}


#include "CalculatorGeneric.h"



void Calculator::Rand(){  // we are not aiming for perfection here.....
	if (stack.size() > 0) {

		int sz =  DIGITS * (uint)( stack.back()->size()-1);
		int s = stack.back()->size()? stack.back()->back() : 0;

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

		BIntPtr temp(new BInt); temp->clear();
		for (; size >= DIGITS; size = size - DIGITS)
			temp->push_back(_Rand(RMOD)); //  _Rand() returns an integer in the range 0..RMOD-1


		int i =  _Rand(RMOD); 

		while (DIGITS > size ) {
			i = i / 10;
			size++;
		}
		i = i% stack.back()->back();

		temp->push_back(i);

		Normalize(temp);

		stack.pop_back();
		stack.push_back(temp);
	}
	else 
		std::cout << "Empty Stack !" << std::endl;

}


void Calculator::Div2(unsigned int Power)
{
	unsigned int _Shift = Power;
	if (stack.size() > 0) {
		int Ctemp = 0;
		int C1 = 0;
		BIntPtr temp(new BInt);
		Pop(*temp);
		int tempSign = temp->size() ? temp->at(0) : 0;
		while (_Shift > 9) // RMOD == 2^9 * 5^9
		{
			Ctemp = 0;
			C1 = 0;
			for (int i = (int) temp->size() - 1; i >= 0; i--)
			{
				Ctemp = temp->at(i) % 512;
				temp->at(i) = temp->at(i) / 512;
				temp->at(i) = (C1 * (RMOD / 512)) + temp->at(i);
				C1 = Ctemp;
			}
			_Shift -= 9;
			if (temp->back() == 0) temp->pop_back();
		}
		while (_Shift > 0)
		{
			Ctemp = 0;
			C1 = 0;
			for (int i = (int) temp->size() - 1; i >= 0; i--)
			{
				Ctemp = temp->at(i) % 2;
				temp->at(i) = temp->at(i) / 2;
				temp->at(i) = (C1 * (RMOD /2)) + temp->at(i);
				C1 = Ctemp;
			}
			_Shift--;
			if (temp->size() && temp->back() == 0) temp->pop_back();
		}
		Normalize(temp);
		if(temp->size()) temp->at(0) |= (tempSign & SIGNMASK) ;
		stack.push_back(temp);
	}
}



std::string* Calculator::ItoA()
{
	char buffer[12];
	std::string* s = new std::string();
	if (stack.size() > 0) {
		BIntPtr o = stack.back(); stack.pop_back();
		for (int i = 0; i < o->size(); i++) {
			sprintf(buffer, FORMATSTRING, (o->at(i)& ~SIGNMASK));
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
		else if (o->at(0) & SIGNMASK)  s->append("-");
		std::reverse(s->begin(), s->end());
		return s;
	}
	s->append("Empty Stack");
	return s;
}



#define OVERALLOCATION 2

void Calculator::Mul()
{

	if (stack.size() < 2)
		std::cout << "Mul() needs two arguments" << std::endl;
	else if (SMALLNUMBERLIMIT && stack[stack.size() - 1]->size() < 1+SMALLNUMBERLIMIT)
		RussianPeasantMult();
	else if (SMALLNUMBERLIMIT && stack[stack.size() - 2]->size() < 1+SMALLNUMBERLIMIT)
		RussianPeasantMult();
	else {
		PrimeFactorDFT pf;
		BIntPtr temp(new BInt); temp->clear();

		int tempSign = (stack[stack.size() - 1]->at(0)  ^ stack[stack.size() - 2]->at(0)) & SIGNMASK;

		factorSeq  factors;

		u64 min_sz = stack[stack.size() - 1]->size() +
			stack[stack.size() - 2]->size();

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
				temp->push_back(t0+t1+t2);//OVERALLOCATION
			}

			Normalize(temp);
			if(temp->size()) temp->at(0) |= tempSign;
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
	for (int ix = 0; ix < A->size(); ix++) {
		int temp = A->at(ix) & ~SIGNMASK;   // temp is  radix 10^9 
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


void Calculator::Div1000000000(BInt& A)
{
	switch (A.size())
	{
	case 0: 
		break;

	case 1: 
		A.pop_back(); 
		break;

	default:
		int sign = A[0] & SIGNMASK;
		for (int ix = 0; A.size() && (ix < (A.size() - 1)); ix++)
		{
			A[ix] = A[ix + 1];
		}
	    A.pop_back();
		A[0] |= sign;
		break;
	}
}

void Calculator::DumpInt(std::string name,  const BInt& arg)
{
	char buffer[12];
	std::string* s = new std::string();
	for (int i = 0; i < arg.size(); i++) {
		sprintf(buffer, FORMATSTRING, arg[i] & ~SIGNMASK);
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
	if (arg[0] & SIGNMASK) s->append("-");
	std::reverse(s->begin(), s->end());
	std::cout << name <<" : " << *s << " ( " << s->length() << " ) " << std::endl;

}

void Calculator::QuotientRemainder()
{
	BIntPtr Quotient(new BInt); Quotient->clear();
	//Quotient->sign = 1;

	int counter = 0;

	if (stack.size() < 2)
		std::cout << "QuotientReminder() needs two arguments" << std::endl;
	else if (IsZero(*stack[stack.size() - 1])) {
		std::cout << "divison by zero" << std::endl;
		return;
	}
	else if ((stack[stack.size() - 1]->size() == 1) &&
		(stack[stack.size() - 2]->size() == 1)) {
		/* small numbers both less than RMOD */
		BIntPtr    Remainder(new BInt); Remainder->clear();
		Quotient->push_back(stack[stack.size() - 2]->at(0) / stack[stack.size() - 1]->at(0));
		Remainder->push_back(stack[stack.size() - 2]->at(0) % stack[stack.size() - 1]->at(0));
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
		int     reciprocal = (RMOD) /(2+ divisor.back());
		BIntPtr _reciprocal(new BInt); _reciprocal->clear(); _reciprocal->push_back(reciprocal);
		int     shift = (int)divisor.size();
		BIntPtr _dividend(new BInt); Dup(*_dividend, dividend);
		BIntPtr _divisor(new BInt);  Dup(*_divisor, divisor);

		stack.push_back(_reciprocal);
		stack.push_back(_dividend);
		Mul(); 	
		if (stack.back()->size())  for (int i = 0; i < shift;i++)
			DIVRMOD(*stack.back());	

		while(1)
		{

			Dup(*Quotient, *stack.back()); 
			stack.push_back(_divisor);
			Mul(); 
			ChangeSign();
			stack.push_back(_dividend); 
			Add();
			if (IsAbiggerNummerically(divisor, *stack.back()))
				break;

			stack.push_back(_reciprocal);
			Mul();
			for (int i = 0; i < shift;i++) {
				DIVRMOD(*stack.back()); 
			}
			if (IsZero(*stack.back())) {
				stack.pop_back();
				Push(1);
			};
			stack.push_back(Quotient);
			Add();
		}
		/* we are done */
		if (stack.back()->size() && (stack.back()->at(0)& SIGNMASK))
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


