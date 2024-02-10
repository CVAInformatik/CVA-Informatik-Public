#pragma once
/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.

*/
typedef unsigned long long uint;

class SHR {

public:

    SHR(uint _taps, uint _init) :taps(_taps), init(_init), state(init)
	{
		/* build mask */
		mask = _taps;
		mask = mask | (mask >> 1);
		mask = mask | (mask >> 2);
		mask = mask | (mask >> 4);
		mask = mask | (mask >> 8);
		mask = mask | (mask >> 16);
		mask = mask | (mask >> 32);

	};

    uint GetState() { return state;	};

	uint GetNextState() { (void) GetNextBit(); return state; };

    bool GetNextBit() {
		uint  t = mask & state & taps;

		unsigned int   w = 0;
		while (t) {
			t = t & (t - 1);
			w++;
		}
		
		if (state == 0)
			state = mask & init;
		else
			state = mask & (state << 1);

		if (w & 1) 	state = state | 1;

		return (w & 1) == 1;		
	};

private:

    uint taps;
    uint mask;
    uint init;
    uint state;
};

/* 
 
 
#include <iostream>
#include "SHR.h"

void testSHR(uint tap)
{
	SHR  shr(tap, 1);
	for (int b = 0; b < 80; b++)
	{
		for (int c = 0; c < 80; c++)
			if (shr.GetNextBit())
				std::cout << "1";
			else
				std::cout << "0";
		std::cout << std::endl;
	}

}


int main()
{
	testSHR(0xE10000);// period 16,777,215

	std::cout << std::endl;
}

*/