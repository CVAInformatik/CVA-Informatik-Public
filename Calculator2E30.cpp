/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/


#include "Calculator2E30.h"
//
#define CALCULATOR Calculator2E30 
#define BINT BInt2E30
#define BINTPTR BInt2E30Ptr
#define DIVRMOD(x)  Div2E30(x)

//#include "CalcUtil.h"

Calculator2E30::Calculator2E30() {
	stack.reserve(10);
	dist = new std::uniform_int_distribution<uint>(0, RMOD - 1);
}

Calculator2E30::~Calculator2E30() {
	delete dist;
}

#include "CalculatorGeneric.h"

// stack operations

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


/* Quick and Dirty and slow, since my previous hack, didn't work*/
std::string* Calculator2E30::ItoA()
{
#define FORMATSTRING2E30 "%09d"

	char buffer[12];
	std::string* s = new std::string();
	int  f = 0;
	if (stack.size() > 0) {
		do {
			Push(1000000000);
			QuotientRemainder();
			f = stack.back()->number[0];
			sprintf(buffer, FORMATSTRING2E30, f);
			char* c1 = buffer;
			char* c2 = buffer + 8;
			while (c1 < c2)
			{
				char t = *c1;
				*c1 = *c2;
				*c2 = t;
				c1++; c2--;
			}

			s->append(buffer);

			Pop();
		} while (!IsEqual(0));

		while (s->size() && (s->back() == '0')) s->pop_back();
		if (s->size() == 0) s->append("0");
		if (stack.back()->sign == -1) s->append("-");
		stack.pop_back();
		std::reverse(s->begin(), s->end());
		return s;
	}
	s->append("Empty Stack");
	return s;

}

