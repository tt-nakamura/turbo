// uses NTL
//   http://www.shoup.net/ntl

#ifndef __ConvolCode_h__
#define __ConvolCode_h__

#include<NTL/GF2X.h>
#include<NTL/mat_GF2.h>
using namespace NTL;

struct ConvolCode {
    static const long coeff[2];
    long m;// memory order
    GF2X g[2];// generator polynomials
    Mat<GF2> X;// table of output bits
    ConvolCode(const long* =coeff);
    void encode(Vec<GF2>&, const Vec<GF2>&);
    long decode(Vec<GF2>&, const Vec<double>&, long=1, double=0);
    void BCJR(Vec<double>& p, const Mat<double>&);
};

void GF2XFromLong(GF2X&, long);
void AddNoise(Vec<double>&, const Vec<GF2>&, double);

#endif // __ConvolCode_h__
