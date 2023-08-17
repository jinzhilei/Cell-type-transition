//
//  Random.h
//  StemCell
//
//  Created by Rongsheng Huang on 2022/12/1.
//

#ifndef Random_h
#define Random_h



#define RANDMAX 23371.0
#define RANDPI 3.1415926

class Random{
private:
    double x[97];
    double c;
public:
    Random()
    {
        int i;
        time_t t;
        double f;
        srand( (unsigned) time(&t));
        for (i = 0; i < 97; i++)
        {
            f = 1.0*rand();
            x[i] = (1.0*fmod(f,RANDMAX))/RANDMAX;
        }
        c = (1.0*fmod(rand(),RANDMAX))/RANDMAX;
    };
    Random(unsigned seed)
    {
        Initialized(seed);
    };
    double Value() // Get the value with uniform distribution;
    {
        double x0,U;
        int r = 0, s = 64;
        int i;
        double d = 7654321.0/16777216.0, d0 = 1677213.0/1677216.0;
        double f;
        if (x[r] >= x[s])
        {
            x0 = x[r] - x[s];
        }
        else
        {
            x0 = x[r] - x[s] + 1;
        }
        if ( c >= d)
        {
            c = c - d;
        }
        else
        {
            c = c - d + d0;
        }
        if (x0 >= c)
        {
            U = x0 - c;
        }
        else
        {
            U = x0 - c + d0;
        }
        for (i = 0; i < 96; i++)
        {
            x[i] = x[i+1];
        }
        x[96] = fmod(x0,1);
        c = fmod(c,1);
        f = fmod(U,1);
        return(f);
    };
    double NormalDistribution()// Get the value with standard normal distribution;
    {
        double f1,f2,x;
        f1 = Value();
        f2 = Value();
        x = sqrt(-2*log(f1))*cos(2*RANDPI*f2);
        return(x);
    };
    double NormalDistribution(double mu, double sigma)// Get the value with normal distribution with mean mu and standard deviation sigma
    {
        double x;
        x = NormalDistribution();
        x = sigma*x+mu;
        return(x);
    };
    int BinomialDistribution(int n, double p)// Get a random number with binormal distribution B(n,p);
    {
        double r;
        double a0, a1, sum;
        int k;
        r = Value();
        if (p > 1.0 - 1e-6)
        {
            k=n;
        }
        else
        {
            a0 = pow(1-p,n);
            sum = a0;
            if (sum > r){
                k=0;
            }
            else{
                for (k=1; k<=n; k++)
                {
                    a1 = a0 * (1.0 * (n-k + 1.0)/(1.0 * k))*(p/(1-p));
                    sum = sum + a1;
                    if (sum > r)
                    {
                        break;
                    }
                    else
                    {
                        a0 = a1;
                    }
                }
            }
        }
        return(k);
    };
    double operator()()// Get a random value;
    {
        return(Value());
    };
    double Value(double a, double b)// Obtain a random value from the interval [a,b]
    {
        return(a + (b-a)*Value());
    };
    double operator()(double a, double b)
    {
        return(Value(a,b));
    };
    
    double BetaDistribution(double a, double b)// Get a random number with BetaDistribution
    {
        double x, y;
        x = GammaDistribution(a,1);
        y = GammaDistribution(b,1);
        return(x/(x+y));
    };
    double GammaDistribution(double a, double b)// Get a random number with GammaDistribution
    {
        int k,n;
        double delta, U, V, W, xi, eta;
        double E=2.71828;
        n=floor(a);
        delta = a - n;
        xi=0;
        if(delta>1e-6)
        {
            for (k=1; k<=100; k++)
            {
                U=Value();
                V=Value();
                W=Value();
                if(U<=E/(E+delta))
                {
                    xi = pow(V,1/delta);
                    eta = W * pow(xi, delta - 1);
                }
                else
                {
                    xi = 1 - log(V);
                    eta = W * pow(E,-xi);
                }
                if (eta < pow(xi, delta-1)*pow(E, -xi))
                {
                    break ;
                }
            }
        }
        for (k=1;k<=n;k++)
        {
            xi = xi - log(Value());
        }
        return(b*xi);
    };
    void Initialized(unsigned seed)// Initialzation of the random number
    {
        int i;
        double f;
        srand(seed);
        for (i = 0; i < 97; i++)
        {
            f = 1.0*rand();
            x[i] = (1.0*fmod(f,RANDMAX))/RANDMAX;
        }
        c = (1.0*fmod(rand(),RANDMAX))/RANDMAX;
    };
    
    double WienerGen()// Generate a Wiener process;
    {
        double f1,f2,x;
        f1 = Value();
        f2 = Value();
        x = sqrt(-2*log(f1))*cos(2*RANDPI*f2);
        return(x);
    };
};

#endif /* Random_h */
