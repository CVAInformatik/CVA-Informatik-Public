
CC = g++
CFLAGS = -g 
# CPPFLAGS = 

%.o  :  %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

PrimeFactorFFT.o : PrimeFactorFFT.cpp

PrimeFactorDFT.o : PrimeFactorDFT.cpp

SlowFFT.o  : SlowFFT.cpp

Calculator.o : Calculator.cpp 

PrimeFactorFFT :  PrimeFactorFFT.o PrimeFactorDFT.o Calculator.o SlowFFT.o
