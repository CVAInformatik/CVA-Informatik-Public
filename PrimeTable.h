#pragma once
/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/


class PrimeTable {

public:
	PrimeTable(unsigned int size) {
		table = new unsigned char[(size / 20) + 1];
		for (unsigned int i = 0; i < (size / 20) + 1; i++) table[i] = 0xFF;
		unsigned int s = 3;
		do {
			if (IsPrime(s)) {
				for (unsigned int i = s + s; i < size ; i = i + s) {
					if ((i & 1) == 0) continue;
					if ((i % 5) == 0) continue;
					unsigned char bm = BitMask(i);
					unsigned int off = i / 20;
					table[off] = table[off] & (~bm);
				}
			}
			s += 2;
		} while (s * s < size);// will terminate after square-root(size) iterations
	
	};

	~PrimeTable() {
		delete[] table;
	};

	bool IsPrime(unsigned int n) {
		if (n <= 1) return false;  // 0,1 are not prime
		if ((n & 1) == 0) return n == 2; // 2 is the oddest prime
		if ((n % 5) == 0) return n == 5; // 5 is prime
		unsigned char bm = BitMask(n);
		unsigned int offset = n / 20;
		return (table[offset] & bm) != 0;
	}

private:

	unsigned char BitMask(unsigned int n){
		switch (n % 20) {
		case  1: return 1; break;
		case  3: return 2; break;
		case  7: return 4; break;
		case  9: return 8; break;
		case 11: return 16; break;
		case 13: return 32; break;
		case 17: return 64; break;
		case 19: return 128; break;
		default: return 0; break; // we should never end up here
		}
	}

	unsigned char *table;
};
