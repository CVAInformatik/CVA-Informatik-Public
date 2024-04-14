// GrayCodeGenerator.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include "GrayCodeGenerator.h"


#define WIDTH 8

void putBit(int w, u64 x)
{
    if (w>1)
    {
        putBit(w - 1, x >> 1);
    }
    printf("%c", x & 1 ? '1' : '0');

}

#ifdef OS_WINDOWS
#define FSTRING "I64"
#else
#define FSTRING "l"
#endif

int main()
{
    GrayCodeGeneratorClass GC;

    GC.Init(WIDTH, 0x7FFFFFFAA);

    u64 o= 0, oold= 0;
    for (u64 i = 0; i <= (((u64) 1) << WIDTH); i++) {
        o = GC.Next();
        printf(" %8" FSTRING "d  %8" FSTRING "X  ",i, o );
        putBit(WIDTH, o);
        printf("  ");
        putBit(WIDTH, o^oold);
        oold = o;
        printf("\n");
    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
