
CC = g++
CFLAGS = -g 
CPPFLAGS =  -O3

%.o  :  %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@


clean:
	rm *.o

CalcUtil.o : CalcUtil.cpp  CalcUtil.h

PrimeFactorFFT.o : PrimeFactorFFT.cpp

PrimeFactorDFT.o : PrimeFactorDFT.cpp

SlowFFT.o  : SlowFFT.cpp

Calculator.o : Calculator.cpp 

PrimeFactorFFT :  PrimeFactorFFT.o PrimeFactorDFT.o Calculator.o SlowFFT.o CalcUtil.o
