
#include "CalculatorType.h"
#include "CalcUtil.h"
#include "Calculator2E30.h"


void ModularExponentiation(BINT& Res, const BINT &a, const BINT &exp, const BINT &mod)
{
	CALCULATOR Square;
    CALCULATOR Prod;
    CALCULATOR Iterator;
    BINT temp;

	Square.Push(a);
    Prod.Push(1);
    Iterator.Push(exp);

	while (!Iterator.IsEqual(0)) {
		if (!Iterator.IsEven()) {
			Square.Dup();
			Square.Pop(temp);
			Prod.Push(temp);
			Prod.Mul(mod);
        }
		Iterator.Div2();
        if (!Iterator.IsEqual(0)) {
            Square.Square(mod);
        }
	}
	Prod.Pop(Res);
}

void ModularMultiplication(BINT& ab, const BINT& a, const BINT& b, const BINT& mod)
{
	CALCULATOR cal;
	cal.Push(a);
	cal.Push(b);
	cal.Mul(mod);
	cal.Pop(ab);
}
void ModularAddition(BINT& aplusb, const BINT& a, const BINT& b, const BINT& mod)
{    
    CALCULATOR cal;

	cal.Push(a);
	cal.Push(b);
	cal.Add(mod);
	cal.Pop(aplusb);
}


void ModularSquare(BINT& Res, const BINT& a, const BINT& mod)
{
    CALCULATOR cal;

    cal.Push(a);
    cal.Square(mod);
    cal.Pop(Res);
}


#define LIMIT 1000000
#define WCOUNT 400

/* we check if M is a prime */
void SquareRootModPrime(BINT& Res, BINT& A, BINT& M)
{
    CALCULATOR cal;
    BINT mtemp;

    PrimeTable pt(LIMIT);

    std::vector<unsigned int>  witnesses;

    witnesses.push_back(2);
    unsigned int i1 = 3;
    while (witnesses.size() < WCOUNT) {
        if (pt.IsPrime(i1)) witnesses.push_back(i1);
        i1 = i1 + 2;
    }

    if (!MillerRabin(M, witnesses)) {
        cal.Push(M);
        std::cout << "probably composite: " << *cal.ItoA() << std::endl;
        return;
    }
    /* we assume M is prime now */
    if ((M.number[0] % 4) == 3) {
        cal.Push(M);
        cal.Push(1);
        cal.Add();
        cal.Div2(2);
        cal.Pop(mtemp);
        ModularExponentiation(Res, A, mtemp, M);
        return;

    }
    SquareRootModM(Res, A, M);
}

/* without checking primality of M*/
void SquareRootModM(BINT& Res, BINT& A, BINT& Mod)
{
    CALCULATOR cal;
    /* */
    BINT Q;
    BINT M;
    BINT P;

    cal.Push(Mod);
    cal.Push(A);
    cal.Jacobi();
    if (!cal.IsEqual(1)) {
        cal.Push(A);
        std::cout << "not a square :"  << *cal.ItoA() << std::endl;
        return;
    }
    cal.Clear();
    /* we assume A is a square and M is a prime
    *  and try to do the Tonneli-Shanks algorithm
    */
    cal.Push(Mod);
    cal.Pop(P);
    cal.Push(P);
    cal.Push(-1);
    cal.Add();
    int S = 0;
    while (!cal.IsEqual(0) && cal.IsEven()) {
        S++;
        cal.Div2();
        if (cal.IsEqual(0)) {
            std::cout << "not a square "  << std::endl;
            return;
        }
    }
    cal.Pop(Q);
    cal.Push(S);
    cal.Pop(M);
    BINT z;
    BINT c;
    BINT t;
    BINT temp, R;
    BINT t1;
    BINT b, temp1;

    do {
        cal.Push(P);
        cal.Dup();
        cal.Rand();
        cal.Pop(z);
        cal.Push(z);
        cal.Jacobi();
    } while (!cal.IsEqual(-1));
    //cal.Push((char*)"20987094688876480237009540038342951804750719680990270636947783175614");
    //cal.Pop(z);
    cal.Clear();
    ModularExponentiation(c, z, Q, P);

    ModularExponentiation(t, A, Q, P);

    cal.Push(Q);
    cal.Push(1);
    cal.Add();
    cal.Div2();
    cal.Pop(temp);
    ModularExponentiation(R, A, temp, P);
loop:
    cal.Push(t);
    if (cal.IsEqual(0)) {
        cal.Pop(Res);
        std::cout << "Root is 0" << std::endl;
        return;
    }
    if (cal.IsEqual(1)) {
        cal.Push(R);
        cal.Pop(Res);
        std::cout << "Root is " << *cal.ItoA() << std::endl;
        return;
    }
    int i = 0;
    do {
        cal.Pop(t1);// t1 = t 
        i++;
        ModularSquare(t1, t1, P);
        cal.Push(t1);
    } while (!cal.IsEqual(1));
    cal.Pop();
    cal.Push(2);
    cal.Push(M);
    cal.Push(i+1);
    cal.ChangeSign();
    cal.Add();
    cal.Exp();
    cal.Pop(temp1);
    ModularExponentiation(b, c, temp1, P);
    cal.Push(i);
    cal.Pop(M);
    ModularSquare(c, b, P);
    ModularMultiplication(t, t, c, P);
    ModularMultiplication(R, R, b, P);
    goto loop;

}

bool MillerRabin(BINT& number, const std::vector<unsigned int>& witnesses)
{
    CALCULATOR cal;

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
                isZ = cal.IsEqual(0);
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
                cal.PushStore("y");
                cal.Square();
                //cal.Dup();
                //std::cout << "y*y " << *cal.ItoA() << std::endl;
                cal.PushStore("m");
                //cal.Dup();
                //std::cout <<"m "<< * cal.ItoA() << std::endl;
                cal.Mod();
                cal.PopStore("y");
                cal.PushStore("y");
                //cal.Push(1);
                if (cal.IsEqual(1)) {
                    cal.PushStore("x");
                    if (!cal.IsEqual(1)) {
                        // x == 1
                        cal.PushStore("m-1");
                        if (!cal.IsEqual()) {
                            //bingo
                            cal.Pop();cal.Pop();//cal.Pop();
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
            //cal.Push(1);
            if (!cal.IsEqual(1)) {
                cal.Pop(); cal.Pop();
                //std::cout << "Composite" << std::endl;
                return false;
            }
        }
        return true;
    }
}

void Factoring(char c[])
{
    CALCULATOR cal;
    BINT temp;

    PrimeTable pt(LIMIT);

    std::vector<unsigned int>  witnesses;

    witnesses.push_back(2);
    unsigned int i = 3;
    while (witnesses.size() < WCOUNT) {
        if (pt.IsPrime(i)) witnesses.push_back(i);
        i = i + 2;
    }

    cal.Push((char*)c);
    cal.Pop(temp);
    if (MillerRabin(temp, witnesses)) {
        cal.Push(temp);
        std::cout << "probably prime: " << *cal.ItoA() << std::endl;
        return;
    }
    // OK we think this is composite
    // let's do some GCD tests
    cal.Push(temp);
    cal.Push(1);
    for (u64 x = 3; x < LIMIT; x = x + 2) {
        if (!cal.IsLarger()) {
            if (pt.IsPrime((unsigned int)x)) {
                cal.Push((int)x);
                cal.Mul();
            }
        }
        else
        {
            cal.Push(temp);
            cal.GCD();
            cal.Push(1);
            if (cal.IsEqual()) {
                // nope not a factor
                cal.Clear();
                cal.Push(temp);
                cal.Push((int)x);
            }
            else {
                cal.Pop();
                std::cout << "factor found " << *cal.ItoA() << std::endl;
                return;
            }
        }
    }

}

void Faculty(BINT& res, int  a)
{
    CALCULATOR cal;

    cal.Push(1);
    for (int i = 2; i <= a; i++)
    {
        cal.Push(i);
        cal.Mul();
    }
    cal.Pop(res);
}

void Convert10E9to2E30(BInt2E30& dest, BInt& src)
{
    Calculator cal;
    BInt  temp;
    BInt  temp1;

    cal.Push(2);
    cal.Push(15);
    cal.Exp();
    cal.PopStore("_2E15");
    cal.Push(src);

    do {
        cal.PushStore("_2E15");
        cal.QuotientRemainder();
        cal.Pop(temp);
        temp1.number.push_back(temp.number[0]);
    } while (!cal.IsEqual(0));

    if ((temp1.number.size() % 2) == 1) temp1.number.push_back(0);

    for (int i = 0; i < temp1.number.size(); i = i + 2)
        dest.number.push_back( temp1.number[i] | (temp1.number[i + 1] << 15));

    dest.sign = src.sign;
}

void Convert2E30to10E9(BInt& dest, BInt2E30&  src)
{
    Calculator cal;
    BInt  temp;
    BInt  temp1;

    /* first convert to radix 2E15, so we can use the calculator */
    for (int i = 0; i < src.number.size(); i = i++)
    {
        temp.number.push_back(src.number[i] & 0x7FFF);
        temp.number.push_back((src.number[i]>>15) & 0x7FFF);
    }
    cal.Push(2);
    cal.Push(15);
    cal.Exp();
    cal.PopStore("_2E15");

    cal.Push(0);
    for (u64 i = temp.number.size() ; i > 0; i--)
    {
        cal.PushStore("_2E15");
        cal.Mul();
        cal.Push(temp.number[i-1]);
        cal.Add();
    }
    cal.Pop(dest);
    dest.sign = src.sign;
}

void MersenneBInt2E20(BInt2E30& dest, uint N)
{
    uint Nt = N;
    dest.number.clear();
    dest.sign = 1;

    while (Nt > 30) {
        dest.number.push_back(0x3FFFFFFF);
        Nt = Nt - 30;
    }

    if (Nt > 0)
    {
        int it = 1;
        Nt--;

        while (Nt-- > 0)  it = it | it << 1;
        dest.number.push_back(it);
    }

}

