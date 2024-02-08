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

typedef unsigned __int64  u64;

typedef  struct _BInt
{
	int sign;
	std::vector<int>  number;
	_BInt() { sign = 1; number.clear(); /*number.push_back(0);*/ };
} BInt;

typedef std::shared_ptr<BInt> BIntPtr;

/* 
		SMALLNUMBERLIMIT :
		0     No small number optimization, FFT all the way down
		1..2  limit is 10^9*SMALLNUMBERLIMIT 
		values beyound 2 will not work 
*/
#define SMALLNUMBERLIMIT 1 

//#define DUMPINT(x,y) DumpInt(x,y)
#define DUMPINT(x,y) 

class Calculator 
{
public:

	Calculator();
	~Calculator();
	// stack operations
	void Push(int i);  // push integer on stack
	void Push(char *c);  // push integer in ascii on stack
	void Push(const BInt &b);  // push BInt (it will be copied)
	inline void Pop() { if (stack.size() > 0) stack.pop_back(); };  //remove TOS 
	inline void Pop(BInt &b) {	if (stack.size() > 0) {	Dup(b, *stack.back());	stack.pop_back();}};  //remove TOS leave a copy in b
	void Swap();  // interchange the two top-most items on stack
	inline void Dup() { if (stack.size() > 0) stack.push_back(stack.back()); };  // push a copy of TOS on Stack
	/* Output     */
	std::string* ItoA(); // pop TOS and return  as a string
	std::string* PrintTOS() { Dup(); return ItoA(); }; // non-destructive print of TOS 
    /* Arithmetic */
	void Mul();  // replace two toplevel items with their product
	void QuotientRemainder();//  replaces the two top elements with Q and R, R < Q
	                        //  such that S1 = Q * S2 + R
	//   
	inline void Square() { Dup(); Mul(); };
	inline void Mod() { QuotientRemainder(); Swap(); Pop();	}
	inline void Div() { QuotientRemainder(); Pop(); }

	void GCD(); // relplace A, B on the stack with [GCD,MA,MB,....] 
	            // such that GCD = A * MA - B * MB
	void ChangeSign(); // replace TOS with -1 * TOS
	void Add();  // replace two toplevel items with their sum
	void Exp(); //  replace   [A,B,....  with [B^A,......
	void Jacobi(); // replace the two toplevel items with the Legendre/Jacobi symbol 
	               //  [A,M,....   -> [(A/M),....  value is -1,0 or 1....
	void Div2(unsigned int Power = 1); // replace the TOS with TOS/2^Power
	void Rand();   // replace the TOS, with a pseudo-Randon number in the closed interval [0...TOS-1]

	/* store operations  */
	inline void PopStore(const std::string& loc) { if (stack.size() > 0) {	Store[loc] = stack.back();	stack.pop_back();}};  //remove TOS and save it at named loc
	void PushStore(const std::string& loc); // push content of named loc on Stack
	void ClearStore(const std::string& loc);// clear named location
	inline void ClearStore() {	Store.clear();	};// clear all locations
	
	/*  Various Predicates */
	int  TOSStatus();//Does not modify the stack return -1,0,1 if the TOS element is respectively negative, zero or postive 
	int  TOSSize();//Does not modify the stack return the size of the TOS BInt 
	bool IsLarger() ;//Does not modify the stack true if TOS is larger than the number below. 
	bool IsEqual();//Does not modify the stack true if TOS is equal to number below. 
	bool IsZero();// Does not modify the stack true if TOS is equal to zero
	bool IsEven();// Does not modify the stack true if TOS is even
	
	/* other */
	inline int  StackSize() { return  (int)stack.size();	}; // how deep is the stack...
	inline int  StoreSize() { return  (int)Store.size(); }; // how many items in the store...
	
	// for internal use......
	void dumpStack(int);
private:
	std::vector<BIntPtr> stack;
	std::map<std::string, BIntPtr> Store;
	bool IsZero(BInt arg);
	bool IsAbiggerNummerically(BIntPtr A, BIntPtr B );
	bool IsAbiggerNummerically(BInt A, BInt B);
	bool IsEqual(BInt A, BInt B);
	bool IsALarger(BInt A, BInt B);

	/* Schönhage-Strassen Helpers */
	void LoadFFT(BIntPtr A, Data* Buffer);
	void Carry(__int64 size, Data* Buffer);

	void GCDAux(BIntPtr A, BIntPtr B);
	void Div2(BInt& A);
	void Div1000000000(BInt& A);
	void Mul2(BInt& A);
	inline void Dup(BInt& D, const BInt& S);
	inline void  Add(BInt& D, const  BInt& S) {BInt temp;	_Add(temp, D, S);Dup(D, temp);	}	;
	inline void _Add(BInt& temp, const BInt& D, const BInt& S) { Push(D);Push(S);Add();	Pop(temp); };
	inline void  Sub(BInt& D,const  BInt& S) { Push(D); Push(S); ChangeSign();	Add();	Pop(D); };

	void Normalize(BInt& b);
	inline void Normalize(BIntPtr b) {	Normalize(*b); };
	void DumpInt(std::string name, const BInt& arg);
	void RussianPeasantMult();
	void RussianPeasantMultAux(int sign, __int64 A, const BInt& B);
	uint  _Rand(uint UpperBound ); // _Rand returns a number in the range 0..UpperBound - 1
	std::random_device rd;
	std::uniform_int_distribution<uint>* dist;
};