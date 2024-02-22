/*
Copyright  © 2024 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
This General Public License does not permit incorporating your program into proprietary programs.If your program is a subroutine library, you may consider it more useful to permit linking proprietary applications with the library.
If this is what you want to do, use the GNU Library General Public License instead of this License.
*/
// for performance measurement
#include "PrimeFactorDFT.h"
#ifdef PERF
#include "windows.h"
#include "profileapi.h"
#endif

#include "SlowFFT.h"
#include <iostream>
#include <random>
#include "Calculator.h"
#include "PrimeTable.h"
#include "CalcUtil.h"

static std::random_device rd;
static std::mt19937 mt(rd());

#define LIMIT 1000000
#define WCOUNT 300


void ClearData(s64 Length, Data* dreal, Data* dimag)
{
    for (s64 i = 0; i < Length; i++) {
        dreal[i] = 0;
        dimag[i] = 0;
    }
}

void test1()
{
    PrimeFactorDFT pf;

    factorSeq  factors;

    std::cout << "Test1 begin " << std::endl;

    //factors.push_back(2);
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);
    //factors.push_back(11);
    //factors.push_back(19);

    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* real = new Data[pf.Status()];
        Data* imag = new Data[pf.Status()];
        for (int i = 0; i < pf.Status(); i++)
        {
            if (i % 4 == 1) real[i] = 1.0;
            else if (i % 2 == 0) real[i] = 0.0;
            else real[i] = -1.0;
            imag[i] = 0.0;
        }

        pf.forwardFFT(real, imag);
        for (int i = 0; i < pf.Status(); i++)
            std::cout << i << " : " << real[i] << "  " << imag[i] << std::endl;

        pf.InverseFFT(real, imag);
        for (int i = 0; i < pf.Status(); i++)
            std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] real;
        delete[] imag;
    }
    std::cout << "Test1 end " << std::endl << std::endl;

}

/* 'balanced' representation */
#define AT 3.0
#define TA (-3.0)
#define GC 1.0
#define CG (-1.0)

void InitDNA(s64 Length, Data *DNAreal, Data *DNAimag)
{
    s64  len = Length / 8;

    for (s64 i = 0; i < Length; i++) { 
        DNAreal[i] = 0; 
        DNAimag[i] = 0;
    }

    for (s64 i = 0; i < len; i++)    {
        DNAreal[i] = AT;
        DNAreal[i + len] = TA;
        DNAreal[i + len + len] = GC;
        DNAreal[i + len + len + len] = CG;
    }

    // DNAreal now contains #len AT, #len TA, #len GC, #len CG entries,
    // followed by slightly more than 4xlen zeroes
    // we permute the first 4 x len entries, to generates some pseudo-DNA

#define PERMCOUNT 10

    static std::uniform_int_distribution<uint>* dist;
    dist = new std::uniform_int_distribution<uint>(0, (4 * ((int)len)) - 1);
    for (uint x = 0; x < PERMCOUNT; x++)
       for (uint ix = 0; ix < 4 * len; ix++)
       {
           uint i1 = 0;
           i1 = dist->operator()(mt);
           Data t = DNAreal[ix];
           DNAreal[ix] = DNAreal[i1];
           DNAreal[i1] = t;
       };

    

}

uint InitSubDNA(s64 substringLength, s64 Length, Data* subDNAreal, Data* subDNAimag, Data *DNAreal)
{

    for (s64 i = 0; i < Length; i++) {
        subDNAreal[i] = 0;
        subDNAimag[i] = 0;
    }

    static std::uniform_int_distribution<uint>* dist;

    dist = new std::uniform_int_distribution<uint>(0, (4 * ((int) (Length/8))) - 1);
    uint i1 = 0;
    do {
        i1 = dist->operator()(mt);
    } while (i1 >= ((4 * ((int)(Length / 8))) - substringLength));
    std::cout << "start:  " << i1 << std::endl;
    /* we swap direction ! */
    s64 iy = Length - 1;
    for (uint ix = 0; ix < substringLength; ix++)
        subDNAreal[iy--] = DNAreal[ix + i1];

    return i1;

}



void Multiply(s64 Length , Data *AxBreal, Data* AxBimag, Data* Areal, Data* Aimag, Data* Breal, Data*  Bimag)
{
    for (s64 i = 0; i < Length; i++) {
        Data tr, ti; 

        tr = Areal[i] * Breal[i] - Aimag[i] * Bimag[i];
        ti = Breal[i] * Aimag[i] + Areal[i] * Bimag[i];
        AxBreal[i] = tr;
        AxBimag[i] = ti;
    }

}

void test2DNA()
{

#define SUBSTRINGLENGTH 300

    PrimeFactorDFT pf;

    factorSeq  factors;


    std::cout << "Test2DNA begin " << std::endl;
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);
    factors.push_back(11);
    factors.push_back(19);
    factors.push_back(31);

    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* DNAreal    = new Data[pf.Status()];
        Data* DNAimag    = new Data[pf.Status()];
        Data* subDNAreal = new Data[pf.Status()];
        Data* subDNAimag = new Data[pf.Status()];
        Data* Matchreal  = new Data[pf.Status()];
        Data* Matchimag  = new Data[pf.Status()];

        InitDNA(pf.Status(),DNAreal, DNAimag);
        InitSubDNA(SUBSTRINGLENGTH, pf.Status(), subDNAreal, subDNAimag, DNAreal );


        pf.forwardFFT(DNAreal, DNAimag);
        pf.forwardFFT(subDNAreal, subDNAimag);

        Multiply(pf.Status(), Matchreal, Matchimag, DNAreal, DNAimag, subDNAreal, subDNAimag);
        pf.ScaledInverseFFT(Matchreal, Matchimag);

        Data max = 0;
        Data max2 = 0;
        s64 maxIndex = 0;
        s64 maxIndex2 = 0;
        for (int i = 0; i < pf.Status(); i++) {
            Data val = (Matchreal[i] * Matchreal[i]) + (Matchimag[i] * Matchimag[i]);
            //std::cout << i << " : " << val << std::endl;
            if (val > max) {
                max2 = max;
                max = val;
                maxIndex2 = maxIndex;
                maxIndex = i;
            }
        }
        std::cout << " maxIndex  : " << maxIndex   << " val : " << max  << std::endl;
        std::cout << " maxIndex2 : " << maxIndex2  << " val : " << max2 << std::endl;
        //pf.InverseFFT(real, imag);
        //for (int i = 0; i < pf.Status(); i++)
        //    std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] DNAreal;
        delete[] DNAimag;
        delete[] subDNAreal;
        delete[] subDNAimag;
        delete[] Matchreal;
        delete[] Matchimag;

    }
    std::cout << "Test2DNA End " << std::endl << std::endl;

}

void test3Convolution()
{


    PrimeFactorDFT pf;

    factorSeq  factors;


    std::cout << "Test3Convolution begin " << std::endl;
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);
    factors.push_back(11);
    //factors.push_back(19);
    //factors.push_back(31);

    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* A1real = new Data[pf.Status()];
        Data* A1imag = new Data[pf.Status()];
        Data* A2real = new Data[pf.Status()];
        Data* A2imag = new Data[pf.Status()];
        Data* Resreal = new Data[pf.Status()];
        Data* Resimag = new Data[pf.Status()];

        ClearData(pf.Status(), A1real, A1imag);
        ClearData(pf.Status(), A2real, A2imag);
        ClearData(pf.Status(), Resreal, Resimag);

        A1real[0] = 1.0/15;
        A1real[1] = 12.0/15;
        A1real[2] = 1.0/15;
        A1real[3] = 0.5/15;
        A1real[pf.Status() - 1] = 0.5/15;
        A2real[500] = 1.0;

        pf.forwardFFT(A1real, A1imag);
        pf.forwardFFT(A2real, A2imag);
#define COUNT 5
        Multiply(pf.Status(), Resreal, Resimag, A1real, A1imag, A2real, A2imag);

        for (int i = 1; i < COUNT; i++)
            Multiply(pf.Status(), Resreal, Resimag, A1real, A1imag, Resreal, Resimag);
        
        pf.ScaledInverseFFT(Resreal, Resimag);

        //pf.InverseFFT(Resreal, Resimag);
        for (int i = 0; i < pf.Status(); i++)
            std::cout << i << " " << Resreal[i] << "  " << Resimag[i] << std::endl;

        delete[] A1real ;
        delete[] A1imag ;
        delete[] A2real ;
        delete[] A2imag ;
        delete[] Resreal;
        delete[] Resimag;

    }
    std::cout << "Test3Convolution End " << std::endl << std::endl;

}


void test4()
{
    PrimeFactorDFT pf;
    SlowFFT sft;

    factorSeq  factors;

    std::cout << "Test1 begin " << std::endl;

    /*factors.push_back(2);
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);*/
    factors.push_back(11);
    factors.push_back(13);
    factors.push_back(19);
    factors.push_back(31);
    pf.SetFactors(factors);
    sft.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* real = new Data[pf.Status()];
        Data* imag = new Data[pf.Status()];
        Data* sreal = new Data[pf.Status()];
        Data* simag = new Data[pf.Status()];
        for (int i = 0; i < pf.Status(); i++)
        {
            if (i % 4 == 1) sreal[i] = real[i] = 1.0;
            else if (i % 2 == 0) sreal[i] = real[i] = 0.0;
            else sreal[i] = real[i] = -1.0;
            simag[i] = imag[i] = 0.0;
        }
        std::cout << "PFA begin" << std::endl;
        pf.forwardFFT(real, imag);
        std::cout << "PFA end" << std::endl;

        std::cout << "SFT begin" << std::endl;
        sft.forwardFFT(sreal, simag);
        std::cout << "SFT end" << std::endl;

        for (int i = 0; i < pf.Status(); i++) {
            std::cout << "PFA :" << i << " : " << real[i] << "  " << imag[i] << std::endl;
            std::cout << "Ref :" << i << " : " << sreal[i] << "  " << simag[i] << std::endl << std::endl;
        }
        //pf.InverseFFT(real, imag);
        //for (int i = 0; i < pf.Status(); i++)
        //    std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] real;
        delete[] imag;
        delete[] sreal;
        delete[] simag;
    }
    std::cout << "Test1 end " << std::endl << std::endl;

}

void test5perf()
{
#ifdef PERF
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
#endif

    PrimeFactorDFT pf;

    factorSeq  factors;

    std::cout << "Test5 begin " << std::endl;

    //factors.push_back(2);
    factors.push_back(3);
    factors.push_back(5);
    factors.push_back(7);
    factors.push_back(11);
    //factors.push_back(13);
    //factors.push_back(17);
    factors.push_back(19);
    factors.push_back(31);
    pf.SetFactors(factors);
    std::cout << "status " << pf.Status() << std::endl;

    if (pf.Status() > 0) {

        Data* real = new Data[pf.Status()];
        Data* imag = new Data[pf.Status()];
        for (int i = 0; i < pf.Status(); i++)
        {
            if (i % 4 == 1)  real[i] = 1.0;
            else if (i % 2 == 0)  real[i] = 0.0;
            else real[i] = -1.0;
            imag[i] = 0.0;
        }
        std::cout << "PFA begin" << std::endl;
#ifdef PERF
       QueryPerformanceFrequency(&Frequency);
        QueryPerformanceCounter(&StartingTime);
#endif
        pf.forwardFFT(real, imag);
#ifdef PERF
        QueryPerformanceCounter(&EndingTime);
        ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
        ElapsedMicroseconds.QuadPart *= 1000000;
        ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
        std::cout << "Length " << pf.Status() << "  Elapsed time (microseconds): " << ElapsedMicroseconds.QuadPart << std::endl;
#endif
        std::cout << "PFA end" << std::endl;


       /* for (int i = 0; i < pf.Status(); i++) {
            std::cout << "PFA :" << i << " : " << real[i] << "  " << imag[i] << std::endl;
        }*/
        //pf.InverseFFT(real, imag);
        //for (int i = 0; i < pf.Status(); i++)
        //    std::cout << i << " " << real[i] << "  " << imag[i] << std::endl;

        delete[] real;
        delete[] imag;
    }
    std::cout << "Test5 end " << std::endl << std::endl;

}

void test6Calc()
{
    Calculator c;
    //char num[] = "261261924691694619461924691649164";
    char num[] = "90000000";
    std::cout << num << std::endl;
    c.Push(num);
    c.Dup();
    std::cout << *c.ItoA() << std::endl;
    //char num1[] = "1000000";
    char num1[] = "723972359729357927536526515164159069889121"
        //"72397235972935792753652651516415906988912444" 
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        //"72397235972935792753652651516415906988912444"
        ;

    std::cout << num1 << std::endl;
    c.Push(num1);
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();
    std::cout << *c.ItoA() << std::endl << std::endl;
    c.Dup();
    c.Mul();
    c.Dup();

    std::string *s = c.ItoA();
    std::cout << "size " << s->length()  << std::endl;
    std::cout << *s << std::endl << std::endl << std::endl;
}



void test7Calc()
{
    Calculator c;
    char num[] = "261261924691694619461924691649164";
    c.Push(num);
    c.PopStore("test");
    std::cout << *c.ItoA() << std::endl;
    c.PushStore("test");
    std::cout << *c.ItoA() << std::endl;
    c.PushStore("test1");
    c.Push(0xFFFFFFFF);
    std::cout << *c.ItoA() << std::endl;
    c.Push(0x7FFFFFFF);
    std::cout << *c.ItoA() << std::endl << std::endl;
}

void test8Calc()
{
    Calculator c;

    char num1[] = "500000000";
    char num2[] = "500000001";
    char num3[] = "3996879";
    char num4[] = "4637923";
    char num3000[] = "3996879000";
    char num4000[] = "4637923000";
    char num[] = "261261924691694619461924691649164";


    c.Push(num1);
    c.Push(num2);
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num1);
    c.ChangeSign();
    c.Push(num2);
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num2);
    c.ChangeSign();
    c.Push(num1);
    c.Add();
    std::cout << *c.ItoA() << std::endl;

    c.Push(num3);
    c.Push(num4);
    c.ChangeSign();
    c.Add();
    std::cout << *c.ItoA() << std::endl;

    c.Push(num4);
    c.Push(num3);
    c.ChangeSign();
    c.Add();
    std::cout << *c.ItoA() << std::endl;

    c.Push(num3000);
    c.Push(num4000);
    c.ChangeSign();
    c.Add();
    std::cout << *c.ItoA() << std::endl;

    c.Push(num4000);
    c.Push(num3000);
    c.ChangeSign();
    c.Add();
    std::cout << *c.ItoA() << std::endl;



    c.Push(num);
    c.Dup();
    c.Dup();
    std::cout << *c.ItoA() << std::endl;
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num);
    c.Dup();
    c.Dup();
    std::cout << *c.ItoA() << std::endl;
    c.ChangeSign();
    c.Dup();
    std::cout << *c.ItoA() << std::endl;
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num);
    c.Push(-1);
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(1);
    c.Push(num);
    c.ChangeSign();
    c.Add();
    std::cout << *c.ItoA() << std::endl ;
    std::cout  << std::endl;
}

void test9Calc()
{
    Calculator c;
    char num1[] = "77777777777777777777777777777777";
    char num2[] = "33333333333333333333333333333333";
    c.Push(num1);
    c.Push(num2);
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num1);
    c.Push(num2);
    c.ChangeSign();
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num2);
    c.Push(num1);
    c.ChangeSign();
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num2);
    c.ChangeSign();
    c.Push(num1);
    c.ChangeSign();
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num2);
    c.ChangeSign();
    c.Push(num1);
    c.Add();
    std::cout << *c.ItoA() << std::endl;
    c.Push(num1);
    c.ChangeSign();
    c.Push(num2);
    c.Add();
    std::cout << *c.ItoA() << std::endl  << std::endl;
}

void test10Calc()
{
    Calculator c;
    c.Push(1);
    c.Push(1);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    c.Push(0); 
    c.Push(1);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    c.Push(2);
    for (int i = 0; i < 20; i++) {
        c.Dup();
        c.Mul();
        c.Dup();
        std::string* s = c.ItoA();
        std::cout << "size " << s->length() << std::endl;
        std::cout << *s << std::endl << std::endl;
    }

}

void test11CalcGCD()
{
    Calculator c;
//    char num[] = "693";//"555555555";
//    char num1[] = "609";//"555";
    char num1[] = "555555555";
    char num[] = "555";
    c.Push(num1);
    c.Push(num);
    c.GCD();
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    char num2[] = "5555";
    char num3[] = "557";
    c.Push(num2);
    c.Push(num3);
    c.GCD();
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    char num4[] = "110893008";
    char num5[] = "7448755608";
    c.Push(num4);
    c.Push(num5);
    c.GCD();
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;
    std::cout << *c.ItoA() << std::endl;

}

void test12CalcSmall()
{
    Calculator c;
    //    char num[] = "693";//"555555555";
    //    char num1[] = "609";//"555";
    char num1[] = "777";
    char num[] = "555";
    c.Push(num1);
    c.Push(num);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    char num01[] = "77777777777777";
    char num0[] = "555";
    c.Push(num01);
    c.Push(num0);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    char num21[] = "777";
    char num2[] = "5555555555555555";
    c.Push(num21);
    c.Push(num2);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    char num31[] = "77777777777777777777777777777777777777";
    char num3[] = "555555555555555555555555555555555555555";
    c.Push(num31);
    c.Push(num3);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;
    char num41[] = "7777777777777777777777777777777777777777777777777777777777777777777777777777";
    char num4[] = "555555555555555555555555555555555555555555555555555555555555555555555555555555";
    c.Push(num41);
    c.Push(num4);
    c.Mul();
    std::cout << *c.ItoA() << std::endl;

}

void test13QuotientReminder()
{
    Calculator c;
    
    char num1[] = "777777777777777";
    char num[] = "5555555555555";

    //c.Push(num1);
    //c.Push(num);
    //c.QuotientRemainder();
    //std::cout << "Reminder: " <<  *c.ItoA() << std::endl;
    //std::cout << "Quotient: " <<  *c.ItoA() << std::endl;

    //char num21[] = "777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777";
    //char num2[] = "5555555555555555555555555555555555555555555555555555555555555555555555555555555";
    //c.Push(num21);
    //c.Push(num2);
    //c.QuotientRemainder();
    //std::cout << "Reminder: " << *c.ItoA() << std::endl;
    //std::cout << "Quotient: " << *c.ItoA() << std::endl;

#define LIMIT1 100000
#define STEP  177395769

    c.Push(num1);
    c.PopStore("num1");
    c.Push(num);
    c.PopStore("num");
    for (int i = 0; i < LIMIT1; i++)
    {
        c.PushStore("num1");
        //c.Dup();
        //std::cout << "num1: " << *c.ItoA() << std::endl;
        c.PushStore("num");
        c.Push(i);
        c.Push(STEP);
        c.Mul();
        c.Add();
        //c.Dup();
        //std::cout << "num:      " << *c.ItoA() << std::endl;
        c.QuotientRemainder();
        c.PopStore("Remainder:");
        c.PopStore("Quotient: ");
        //c.PushStore("Remainder:");
        //std::cout << "Reminder: " << *c.ItoA() << std::endl;
        //c.PushStore("Quotient: ");
        //std::cout << "Quotient: " << *c.ItoA() << std::endl;
        c.PushStore("num");
        c.Push(i);
        c.Push(STEP);
        c.Mul();
        c.Add();
        c.PushStore("Quotient: ");
        std::cout << "Quotient: " << *c.PrintTOS() << std::endl;
        c.Mul();
        c.PushStore("Remainder:");
        c.Add();
        c.PushStore("num1");
        if (!c.IsEqual()) {
            c.Swap();
            std::cout << "num * Quotient + Remainder: " << *c.ItoA() << std::endl;
            std::cout << "Num1                      : " << *c.ItoA() << std::endl;
        }

    }
    c.ClearStore();
}

void test14Exp()
{
    Calculator c;
    c.Push(2);
    c.Push(86243);
    c.Exp();
    c.Push(-1);
    c.Add();
    std::cout << "Exp " << *c.ItoA() << std::endl;
    c.Push(2);
    c.Push(13466917);
    c.Exp();
    c.Push(-1);
    c.Add();
    std::cout << "Exp " << *c.ItoA() << std::endl;
    c.Push(2);
    c.Push(20996011);
    c.Exp();
    c.Push(-1);
    c.Add();
    std::cout << "Exp " << *c.ItoA() << std::endl << std::endl;

}

void test15Jacobi()
{
	Calculator c;
	std::cout << "Jacobi ( ";
	c.Push(9907);
	c.Dup();
	std::cout << *c.ItoA() << "/ ";
	c.Push(1001);
	c.Dup();
	std::cout << *c.ItoA() << " ) = ";
	c.Jacobi();
	std::cout << *c.ItoA() << std::endl;


	for (int j = 1; j < 60; j += 2)
		for (int i = 1; i < 31; i++) {
			std::cout << "Jacobi ( ";
			c.Push(j);
			c.Dup();
			std::cout << *c.ItoA() << "/ ";
			c.Push(i);
			c.Dup();
			std::cout << *c.ItoA() << " ) = ";
			c.Jacobi();
			std::cout << *c.ItoA() << std::endl;

		}

}


bool MillerRabin(BInt& number);
bool MillerRabin(BInt& number, const std::vector<unsigned int>& witnesses);


void test16MillerRabin(char c[]) {

    PrimeTable pt(LIMIT);

    std::vector<unsigned int>  witnesses;

    witnesses.push_back(2);
    unsigned int i = 3;
    while (witnesses.size() < WCOUNT) {
        if (pt.IsPrime(i)) witnesses.push_back(i);
        i = i + 2;
    }

    Calculator cal; 
    cal.Push((char*) c );
    for (int itx = 0; itx < 2000; itx++) {
        cal.Push((char*)c);
        cal.Push(2);
        cal.Push(itx);
        cal.Mul();
        cal.Add();    
        BInt temp;
        cal.Pop(temp);
        if (MillerRabin(temp, witnesses)) {
            cal.Push(temp);
            std::cout << "probably prime: " << *cal.ItoA() << std::endl;
        }
    }
}

void test17()
{
    Calculator c;

    c.Push((char *) "9365865165198658618561865816581658");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }
    c.Push((char*)"999999999");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }
    c.Push((char*)"999");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }

    c.Push((char*)"99");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }

    c.Push((char*)"9");
    for (int i = 0; i < 1000; i++) {
        c.Dup();
        c.Dup();
        std::cout << "Arg    " << *c.ItoA() << std::endl;
        c.Rand();
        std::cout << "Rand() " << *c.ItoA() << std::endl;

    }
    c.Push((char*)"0");
    c.Dup();
    c.Dup();
    std::cout << "Arg    " << *c.ItoA() << std::endl;
    c.Rand();
    std::cout << "Rand() " << *c.ItoA() << std::endl;

}

void test18()
{
    Calculator c;

    c.Push((char*)"9999999999999999999999999999999999");
    c.Dup();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.Dup();
    c.ChangeSign();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.Dup();
    c.ChangeSign();
    c.Swap();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.ChangeSign();
    c.Dup();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.Push((char*)"999999");
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.ChangeSign();
    c.Push((char*)"999999");
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.Push((char*)"999999");
    c.ChangeSign();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"9999999999999999999999999999999999");
    c.ChangeSign();
    c.Push((char*)"999999");
    c.ChangeSign();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"99999999");
    c.Dup();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"99999999");
    c.Dup();
    c.ChangeSign();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"99999999");
    c.Dup();
    c.ChangeSign();
    c.Swap();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;
    c.Push((char*)"99999999");
    c.ChangeSign();
    c.Dup();
    c.Mul();
    std::cout << "res    " << *c.ItoA() << std::endl;

}

bool MillerRabin(BInt& number)
{
	PrimeTable pt(LIMIT);

	std::vector<unsigned int>  witnesses;

	witnesses.push_back(2);
	unsigned int i = 3;
	while (witnesses.size() < WCOUNT) {
		if (pt.IsPrime(i)) witnesses.push_back(i);
		i = i + 2;
	}

    Calculator cal;
	cal.Push(number);
	if (cal.IsEven()) {
		std::cout << "argument must be odd " << std::endl;
        return false;
	}
	else {
		cal.PopStore("m"); //store argument in "m"
		cal.PushStore("m");
		cal.Push(-1);
		cal.Add();
		cal.PopStore("m-1");
		cal.PushStore("m-1");
		int s = 0;
		while (cal.IsEven()) {
			cal.Div2();
			cal.PopStore("d");
			cal.PushStore("d");
			s++;
		}
		cal.Pop(); //remove d from stack, it is saved in the store
		// m-1 has now been cleaned of factor 2^s
		for (size_t ix = 0; ix < witnesses.size(); ix++)
		{
			// calculate  x^d mod n 
			cal.Push(witnesses[ix]); // x on stack
			cal.PopStore("Multiplier");// x in multiplier
			bool isZ = false;
			cal.Push(1);
			cal.PopStore("x");
			do {
				cal.PushStore("d");
				bool isE = cal.IsEven();
				cal.Div2();
				isZ = cal.IsZero();
				cal.PopStore("d");
				cal.PushStore("x");
				if (!isE) {
					cal.PushStore("Multiplier");
					cal.Mul();
					cal.PushStore("m");
					cal.Mod();
				}
				cal.PopStore("x");
				if (!isZ) {
					cal.PushStore("Multiplier");
					cal.Square();
					cal.PushStore("m");
					cal.Mod();
					cal.PopStore("Multiplier"); // 
				}
			} while (!isZ);
			// repeat s times...
			cal.PushStore("x");
			cal.PopStore("y");
			for (int sx = 0; sx < s; sx++) {
				//cal.dumpStack(4);
				cal.PushStore("y");
				cal.Square();
				cal.PushStore("m");
				cal.Mod();
				cal.PopStore("y");
				cal.PushStore("y");
				cal.Push(1);
				if (cal.IsEqual()) {
					cal.PushStore("x");
					if (!cal.IsEqual()) {
						// x == 1
						cal.PushStore("m-1");
						if (!cal.IsEqual()) {
							//bingo
							cal.Pop();cal.Pop();cal.Pop();
							//std::cout << "Composite" << std::endl;
                            return false;
						}
						cal.Pop(); //get rid "m-1" on stack
					}
					else
						cal.Pop(); // get rid of "x" on stack
					cal.Pop(); // get rid of "x" on stack
				}
				cal.Pop(); // get rid of 1 on stack
				cal.PopStore("x"); //save x = y

			}
			cal.PushStore("y");
			cal.Push(1);
			if (!cal.IsEqual()) {
				cal.Pop(); cal.Pop();
				//std::cout << "Composite" << std::endl;
				return false;
			}
		}
		return true;
	}

}

void test19()
{
    Calculator cal;
    cal.Push(859433);
    cal.Push(3021377);
    cal.Mul();
    std::string *s = cal.ItoA();
    Factoring((char *) s->c_str());

    Factoring((char*)"19777122841");

}


void test20() 
{
    Calculator cal;
    BInt P;
    BInt A;
    BInt Res;
#define P224 1
#if P224
    // NIST P-224 
    cal.Push((char*)"26959946667150639794667015087019630673557916260026308143510066298881");
    cal.Pop(P);
    cal.Push(2021);
    cal.Dup();
    cal.Dup();
    cal.Mul();
    cal.Mul();
    cal.Push(2021);
    cal.Push(-3);
    cal.Mul();
    cal.Add();
    cal.Push((char*)"18958286285566608000408668544493926415504680968679321075787234672564");
    cal.Add();
    cal.Pop(A);
    SquareRootModM(Res, A, P);
    cal.Push(A);
    std::cout << " A:   " << *cal.ItoA() << std::endl;
    cal.Push(P);
    std::cout << " P:   " << *cal.ItoA() << std::endl;
    cal.Push(Res);
    cal.Dup();
    std::cout << " Res: " << *cal.ItoA() << std::endl;
    cal.Square();
    cal.Dup();
    std::cout << " Res * Res  : " << *cal.ItoA() << std::endl;
    cal.Push(P);
    cal.Mod();
    std::cout << " Res * Res  mod P: " << *cal.ItoA() << std::endl;
    return;

#else 
    cal.Push((char* )"2147483647");
    //cal.Push((char* )"43");
    cal.Pop(P);
    cal.Push((char* ) "3497491");
    //cal.Push((char* ) "6");
    cal.Pop(A);
    cal.Push(P);
    cal.Push(A);
    cal.Jacobi();
    if (cal.IsOne()) {
        SquareRootModPrime(Res, A, P); 
        cal.Push(A);
        std::cout << " A:   " << *cal.ItoA() << std::endl;
        cal.Push(P);
        std::cout << " P:   " << *cal.ItoA() << std::endl;
        cal.Push(Res);
        cal.Dup();
        std::cout << " Res: " << *cal.ItoA() << std::endl;
        cal.Square();
        cal.Dup();
        std::cout << " Res * Res  : " << *cal.ItoA() << std::endl;
        cal.Push(P);
        cal.Mod();
        std::cout << " Res * Res  mod P: " << *cal.ItoA() << std::endl;
        return;
    }
#endif
}



void test21()
{
#ifdef PERF
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;
#endif

    Calculator cal;
    BInt res;
    int arg = 200;

#ifdef PERF
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
#endif

    Faculty(res, arg);
#ifdef PERF
    QueryPerformanceCounter(&EndingTime);
    ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMicroseconds.QuadPart *= 1000000;
    ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
    std::cout <<  "Elapsed time(microseconds) : " << ElapsedMicroseconds.QuadPart << std::endl;
#endif
    cal.Push(res);
    std::cout << " " << arg << "! : " << *cal.ItoA() << std::endl;
}

int main()
{
    //test1();
    //test2DNA();
    //test3Convolution();
    //test4();
    //test5perf();
    //test6Calc();
    //test7Calc();
    //test8Calc();
    //test9Calc();
    //test10Calc();
    //test11CalcGCD();
    //test12CalcSmall();
    //test13QuotientReminder();
    //test14Exp();
    //test15Jacobi();
    //test17();
    //test16MillerRabin((char*)"2147483647");// Mersenne Prime
    //test16MillerRabin((char *) "1228467");
    //test16MillerRabin((char*)"333228469");
    //test16MillerRabin((char*)"19777122847");
    //test18();
    //Factoring((char*)"2147483649");
    //test19();
    //test20();
    test21();
    std::cout << "Done !\n";
}
