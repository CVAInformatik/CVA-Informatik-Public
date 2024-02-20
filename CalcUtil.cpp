

#include "CalcUtil.h"

void ModularExponentiation(BInt& Res, const BInt &a, const BInt &exp, const BInt &mod)
{
	Calculator Square;
	Calculator Prod;
	Calculator Iterator;

	Square.Push(a);
    Prod.Push(1);
    Iterator.Push(exp);
    BInt temp;

	while (!Iterator.IsZero()) {
		if (!Iterator.IsEven()) {
			Square.Dup();
			Square.Pop(temp);
			Prod.Push(temp);
			Prod.Mul(mod);
   //         Prod.Push(mod);
			//Prod.Mod();
        }
		Iterator.Div2();
        if (!Iterator.IsZero()) {
            Square.Square(mod);
            //Square.Push(mod);
            //Square.Mod();
        }
	}
	Prod.Pop(Res);
}

void ModularMultiplication(BInt& ab, const BInt &a, const BInt &b, const BInt &mod)
{
	Calculator cal;
	cal.Push(a);
	cal.Push(b);
	cal.Mul(mod);
    //cal.Push(mod);
    //cal.Mod();
	cal.Pop(ab);
}

void ModularAddition(BInt& aplusb, const BInt &a, const BInt &b, const BInt &mod)
{
	Calculator cal;
//	cal.Push(mod);
	cal.Push(a);
	cal.Push(b);
	cal.Add(mod);
    //cal.Push(mod);
    //cal.Mod();
	cal.Pop(aplusb);
}


void ModularSquare(BInt& Res, const BInt &a, const BInt &mod)
{
    Calculator cal;
    cal.Push(a);
    cal.Square(mod);
    //cal.Push(mod);
    //cal.Mod();
    cal.Pop(Res);
}


#define LIMIT 1000000
#define WCOUNT 400

/* we check if M is a prime */
void SquareRootModPrime(BInt& Res, BInt& A, BInt& M)
{
    PrimeTable pt(LIMIT);

    std::vector<unsigned int>  witnesses;

    witnesses.push_back(2);
    unsigned int i1 = 3;
    while (witnesses.size() < WCOUNT) {
        if (pt.IsPrime(i1)) witnesses.push_back(i1);
        i1 = i1 + 2;
    }

    Calculator cal;
    if (!MillerRabin(M, witnesses)) {
        cal.Push(M);
        std::cout << "probably composite: " << *cal.ItoA() << std::endl;
        return;
    }
    /* we assume M is prime now */
    if ((M.number[0] % 4) == 3) {
        BInt mtemp;
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
void SquareRootModM(BInt& res, BInt& A, BInt& Mod)
{
    Calculator cal;
    /* */
    BInt Q;
    BInt M;
    BInt P;
    cal.Push(Mod);
    cal.Push(A);
    cal.Jacobi();
    if (!cal.IsOne()) {
        cal.Push(A);
        std::cout << "not a square :" << cal.ItoA() << std::endl;
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
    while (!cal.IsZero() && cal.IsEven()) {
        S++;
        cal.Div2();
        if (cal.IsZero()) {
            std::cout << "not a square :" << *cal.ItoA() << std::endl;
            return;
        }
    }
    cal.Pop(Q);
    cal.Push(S);
    cal.Pop(M);
    BInt z;
    BInt c;
    BInt t;
    do {
        cal.Push(P);
        cal.Dup();
        cal.Rand();
        cal.Pop(z);
        cal.Push(z);
        cal.Jacobi();
    } while (!cal.IsMinusOne());
    cal.Clear();
    ModularExponentiation(c, z, Q, P);
    ModularExponentiation(t, A, Q, P);
    cal.Push(Q);
    cal.Push(1);
    cal.Add();
    cal.Div2();
    BInt temp,R;
    cal.Pop(temp);
    ModularExponentiation(R, A, temp, P);
loop:
    cal.Push(t);
    if (cal.IsZero()) {
        cal.Pop(res);
        //std::cout << "Root is 0" << std::endl;
        return;
    }
    if (cal.IsOne()) {
        cal.Push(R);
        cal.Pop(res);
        //std::cout << "Root is " << *cal.ItoA() << std::endl;
        return;
    }
    int i = 0;
    cal.Push(t); 
    BInt t1;
    cal.Pop(t1);// t1 = t 
    do {
        i++;
        ModularSquare(t1, t1, P);
        cal.Push(t1);
    } while (!cal.IsOne());
    cal.Push(2);
    cal.Push(M);
    cal.Push(i+1);
    cal.ChangeSign();
    cal.Add();
    cal.Exp();
    BInt b, temp1;
    cal.Pop(temp1);
    ModularExponentiation(b, c, temp1, P);
    cal.Push(i);
    cal.Pop(M);
    ModularSquare(c, b, P);
    ModularMultiplication(t, t, c, P);
    ModularMultiplication(R, R, b, P);
    goto loop;

}


bool MillerRabin(BInt& number, const std::vector<unsigned int>& witnesses)
{
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
            //std::cout << "Witness: " << witnesses[ix] << std::endl;
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


void Factoring(char c[])
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
    cal.Push((char*)c);
    BInt temp;
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
