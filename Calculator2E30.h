#pragma once
/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/


#include <random>
#include <vector>
#include <string>
#include <vector>
#include  <memory>
#include <algorithm>
#include <iostream>
#include <map>
#include "PrimeFactorDFT.h"

typedef  struct _BInt2E30
{
	int sign;
	std::vector<int>  number;
	_BInt2E30() { sign = 1; number.clear(); };
} BInt2E30;

typedef std::shared_ptr<BInt2E30> BInt2E30Ptr;

#define SMALLNUMBERLIMIT 1 
//#define DUMPINT(x,y) DumpInt(x,y)
#define DUMPINT(x,y) 

//#include "CalcUtil.h"



//  RMOD == 2^30
#define RMOD     0x40000000
#define RMODMASK 0xC0000000
class Calculator2E30
{
public:

	Calculator2E30();
	~Calculator2E30();
	// stack operations
	void Push(int i);  // push integer on stack
	void Push(const BInt2E30& b);  // push BInt (it will be copied)
	inline void Pop() { if (stack.size() > 0) stack.pop_back(); };  //remove TOS 
	inline void Pop(BInt2E30& b) { if (stack.size() > 0) { Dup(b, *stack.back());	stack.pop_back(); } };  //remove TOS leave a copy in b
	void Swap();  // interchange the two top-most items on stack
	void Dup();  // push a copy of TOS on Stack
	inline void Clear() { stack.clear(); }; // if you want a fresh stack. The Store is not changed

	/* ASCII Conversions    */
	void Push(char* c);  // push integer in ascii on stack
	std::string* ItoA(); // pop TOS and return  as a string
	std::string* PrintTOS() { Dup(); return ItoA(); }; // non-destructive print of TOS 

	/* Arithmetic */
	void Mul();  // replace two toplevel items with their product
	void QuotientRemainder();//  replaces the two top elements with Q and R, R < Q
	//  such that S2 = Q * S1 + R
	// 
	void GCD(); // relplace A, B on the stack with [GCD,MA,MB,....] 
	// such that GCD = A * MA - B * MB
	void ChangeSign(); // replace TOS with -1 * TOS
	void Add();  // replace two toplevel items with their sum
	void Exp(); //  replace   [A,B,....  with [B^A,......
	void Jacobi(); // replace the two toplevel items with the Legendre/Jacobi symbol 
	//  [A,M,....   -> [(A/M),....  value is -1,0 or 1....
	void Div2(unsigned int Power = 1); // replace the TOS with TOS/2^Power
	void Rand();   // replace the TOS, with a pseudo-Randon number in the closed interval [0...TOS-1]


	//   Short Cuts, things we do a lot...
	//   Exp(BInt mod) is probably a bad idea, 
	//   use CalcUtil.h's  void ModularExponentiation(BInt& Res, const BInt &a, const BInt &exp, const BInt &mod);
	//   instead
	inline void Square() { Dup(); Mul(); };
	inline void Square(BInt2E30 Mod) { Dup(); Mul(Mod); };
	inline void Mod() { QuotientRemainder(); Swap(); Pop(); }
	inline void Div() { QuotientRemainder(); Pop(); }
	inline void Mul(BInt2E30 Mod) { Mul(); Push(Mod); this->Mod(); }
	inline void Add(BInt2E30 Mod) { Add(); Push(Mod); this->Mod(); }

	// 
	/* store operations  */
	/*inline*/ void PopStore(const std::string& loc);// { if (stack.size() > 0) { Store[loc] = stack.back();	stack.pop_back(); } };  //remove TOS and save it at named loc
	void PushStore(const std::string& loc); // push content of named loc on Stack
	void ClearStore(const std::string& loc);// clear named location
	inline void ClearStore() { Store.clear(); };// clear all locations

	/*  Various Predicates */
	//int  TOSStatus();//Does not modify the stack return -1,0,1 if the TOS element is respectively negative, zero or postive 
	int  TOSSize();//Does not modify the stack return the size of the TOS BInt 
	bool IsLarger();//Does not modify the stack true if TOS is larger than the number below. 
	bool IsEqual();//Does not modify the stack true if TOS is equal to number below. 
	bool IsEqual(int n);// oes not modify the stack true if TOS is equal to number below. 
	bool IsEven();// Does not modify the stack true if TOS is even

	/* other */
	int  StackSize(); // how deep is the stack...
	int  StoreSize(); // how many items in the store...

	// for internal use......
	void dumpStack(int p);

private:
	std::vector<BInt2E30Ptr> stack;
	std::map<std::string, BInt2E30Ptr> Store;
	bool IsZero(const BInt2E30&  arg);
	bool IsAbiggerNummerically(BInt2E30Ptr A, BInt2E30Ptr B);
	bool IsAbiggerNummerically(BInt2E30 A, BInt2E30 B);
	bool IsEqual(BInt2E30 A, BInt2E30 B);
	bool IsALarger(BInt2E30 A, BInt2E30 B);

	void AddAux(int Asign, BInt2E30& A, int Bsign, const BInt2E30& B);

	/* Schönhage-Strassen Helpers */
	void LoadFFT(BInt2E30Ptr A, Data* Buffer);
	void Carry(s64 size, Data* Buffer);

	void GCDAux(BInt2E30Ptr A, BInt2E30Ptr B);
	void Div2(BInt2E30& A);
	void Div2E30(BInt2E30& A);
	void Mul2(BInt2E30& A);
	void Dup(BInt2E30& D, const BInt2E30& S);
	inline void  Add(BInt2E30& D, const  BInt2E30& S) { BInt2E30 temp;	_Add(temp, D, S);Dup(D, temp); };
	inline void _Add(BInt2E30& temp, const BInt2E30& D, const BInt2E30& S) { Push(D);Push(S);Add();	Pop(temp); };
	inline void  Sub(BInt2E30& D, const  BInt2E30& S) { Push(D); Push(S); ChangeSign();	Add();	Pop(D); };

	void Normalize(BInt2E30& b);
	inline void Normalize(BInt2E30Ptr b) { Normalize(*b); };
	void DumpInt(std::string name, const BInt2E30& arg);
	void RussianPeasantMult();
	//void RussianPeasantMultAux(int sign, s64 A, const BInt& B);
	//void SimpleAdditionSubtractionLadder(int sign, s64 A, const BInt& B);
	void SimpleAdditionSubtractionLadder1(int sign, s64 A, int Bsign, const BInt2E30& B);
	uint  _Rand(uint UpperBound); // _Rand returns a number in the range 0..UpperBound - 1
	std::random_device rd;
	std::uniform_int_distribution<uint>* dist;
};