// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/vec_GF2.h>
using namespace NTL;

#define frand() (RandomWord()/ldexp(1, NTL_BITS_PER_LONG))

double gasdev()
// gaussian random deviates of mean 0 and variance 1
{
    static long iset(0);
    static double gset;
    double fac,rsq,v1,v2;

    if(iset) {
        iset = 0;
        return gset;
    }
    do {
        v1 = 2*frand() - 1;
        v2 = 2*frand() - 1;
        rsq = v1*v1 + v2*v2;
    } while(rsq >= 1 || rsq == 0);
    fac = sqrt(-2*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
}

void AddNoise(Vec<double>& d, const Vec<GF2>& c, double sigma)
// add gaussian noise to code c
// sigma = standard deviation of Gaussian noise
// d = corrupted code (probability that c==1), 0<=d<=1
{
    long i;
    double s(sigma*sigma/2);
    d.SetLength(c.length());
    if(s==0) for(i=0; i<c.length(); i++)
        d[i] = IsOne(c[i]);
    else for(i=0; i<c.length(); i++) {
        d[i] = (IsOne(c[i]) ? -1:1) + sigma*gasdev();
        d[i] = 1/(exp(d[i]/s) + 1);
    }
}

void AddNoise(Vec<GF2>& d, const Vec<GF2>& c, double sigma)
// &d==&c is allowed
{
    long i;
    Vec<double> v;
    AddNoise(v,c,sigma);
    d.SetLength(c.length());
    clear(d);
    for(i=0; i<c.length(); i++)
        if(v[i]>=0.5) set(d[i]);// quantize corruption
}