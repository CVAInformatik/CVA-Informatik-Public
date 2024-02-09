PrimeFactorDFT.cpp/h is an implementation of a In-Place-In-Order FFT for length, which are divisible by one or more of these factors: 2,3,5,7,11,13,17,19 and 31.

Calculator.cpp/h is a small Large-integer calculator , which uses signed-magnitude Radix 10^9 numbers to represent large integers. Originally it was intended to demo Schönhage-Strassen multiplikation
and used signed magnitide  radix 10^3 representation, but it turned out to be such fun to use, so I changed it to Radix 10^9 for better space efficieny (the FFT still uses radix 10^3). The Radix  10^9 format
turns out to be pretty convenient for Signed-Magnitude arithmetic, thanks to the 2-bit headroom.

Primetable.h contains a prime-table class, which is an implementation of the well-known compressed prime-table trick and Eratosthene's sieve.

SlowFFT.cpp/h are for validation of the FFT, but not really usable for anything else.

I have tried to unlearn old-style C and adopt a  more modern C++ style, using std-components.

PrimeFactorFFT.cpp contains main() and a number of demo/test routines both for the FFT and the Calculator.

The code is all C++, no assembly required.
