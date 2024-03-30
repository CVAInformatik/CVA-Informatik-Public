#pragma once

/* define if you want Radix 10^9 undefine if you want Radix 2^30 */
//#define CAL10

// for performance measurements on Windows10
//#define PERF



#ifdef CAL10
#define CALCULATOR Calculator 
#define BINT BInt
#define BINTPTR BIntPtr
#define DIVRMOD(x) Div1000000000(x)
#else
#define CALCULATOR Calculator2E30 
#define BINT BInt2E30
#define BINTTPR BInt2E30Ptr
#define DIVRMOD(x)  Div2E30(x)
#endif
