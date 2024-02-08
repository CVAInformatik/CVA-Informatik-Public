/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
This General Public License does not permit incorporating your program into proprietary programs. If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library. 
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/
#include <math.h>

#include "SlowFFT.h"

SlowFFT::SlowFFT(){}

void SlowFFT::SetFactors(factorSeq& _factors){

	state = 1;
	for (int i = 0; i < _factors.size(); i++) state *= _factors[i];

	if (state > 0)
	{
		ResReal  = new Data[state];
		ResImag  = new Data[state];
		CoefReal = new Data[state];
		CoefImag = new Data[state];
	}

	const double PI = 3.141592653589793238463;
	for (int index = 0; index < state; index++)
	{
		CoefReal[index] = cos((index * -2.0 * PI) / state);
		CoefImag[index] = sin((index * -2.0 * PI) / state);
	}
}


void SlowFFT::forwardFFT(Data* real, Data* imag){

	for (int i = 0; i < state; i++)
		ResReal[i] = ResImag[i] = 0.0;

	for (int j = 0; j < state; j++)
	{
		for (__int64 k = 0; k < state; k++)
		{
			__int64 index = j*k;
			index = index % state;

			ResReal[j] += real[k] * CoefReal[index]
				- imag[k] * CoefImag[index];
			ResImag[j] += real[k] * CoefImag[index]
				+ imag[k] * CoefReal[index];
		}
	}

	for (int i = 0; i < state; i++) {
		real[i] = ResReal[i];
		imag[i] = ResImag[i];
	}
}



void SlowFFT::InverseFFT(Data* real, Data* imag){
	forwardFFT(imag, real);
}

void SlowFFT::ScaledInverseFFT(Data* real, Data* imag){

	forwardFFT(imag, real);

	for (int i = 0; i < state; i++) {
		real[i] = real[i] / state;
	    imag[i] = imag[i]/  state ;
	}
}



SlowFFT::~SlowFFT(){

	delete[] ResReal;
	delete[] ResImag;
	delete[] CoefReal;
	delete[] CoefImag;

}
