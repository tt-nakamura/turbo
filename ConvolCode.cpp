// uses NTL
//   http://www.shoup.net/ntl

#include "ConvolCode.h"

#define CONVOL_DEGMAX 8

const long ConvolCode::coeff[2] = {037, 021};

ConvolCode::ConvolCode(const long *c)
// c = array of coefficients of generator polynomials
{
    long i,j;
    for(i=0; i<2; i++)
        GF2XFromLong(g[i], c[i]);
    m = max(deg(g[0]), deg(g[1]));
    if(m > CONVOL_DEGMAX) Error("degree too large");
    X.SetDims(2, 1<<(m+1));
    for(i=0; i<2; i++)
        for(j=1; j<X.NumCols(); j++)
            if(weight(j&c[i])&1) set(X[i][j]);
}

void ConvolCode::encode(Vec<GF2>& t, const Vec<GF2>& s)
// output: t = encoding of s
// input:  s = source data to encode
{
    long i,j,k;
    GF2X f,h;
    conv(f,s);
    t.SetLength((s.length()+m)*2);
    clear(t);
    for(i=0; i<2; i++) {
        mul(h,f,g[i]);
        for(j=0, k=i; j<=deg(h); j++, k+=2)
            t[k] = h[j];// interleave
    }
}

void ConvolCode::BCJR(Vec<double>& p, const Mat<double>& P)
// p[i] = probability that message[i]==1
// input:  p = current guess for probability
//         P = table of transition probability
// output: p = updated probability by BCJR method
{
    long i,j,k;
    long L(p.length()), M(1<<m), N(M<<1);
    double a,b,c;
    Mat<double> A,B,C(P);
    A.SetDims(L,M);
    B.SetDims(L,M);
    for(i=0; i<L; i++)
        for(j=0; j<N; j++)
            C[i][j] *= (j&1 ? p[i] : 1-p[i]);

    A[0][0] = B[L-1][0] = 1;
    for(j=1; j<M; j++) A[0][j] = B[L-1][j] = 0;
    for(i=0; i<L-1; i++) {//forward pass
        for(a=j=0; j<M; j++) {
            k = j|M;
            c = A[i][j>>1]*C[i][j] + A[i][k>>1]*C[i][k];
            A[i+1][j] = c;
            a += c;
        }
        if(a) for(j=0; j<M; j++) A[i+1][j] /= a;// normalize
        else  for(j=0; j<M; j++) A[i+1][j] = 1;
    }
    for(i=L-1; i>0; i--) {//backward pass
        for(b=j=0; j<N; j+=2) {
            k = j&(M-1);
            c = B[i][k]*C[i][j] + B[i][k+1]*C[i][j+1];
            B[i-1][j>>1] = c;
            b += c;
        }
        if(b) for(j=0; j<M; j++) B[i-1][j] /= b;// normalize
        else  for(j=0; j<M; j++) B[i-1][j] = 1;
    }
    for(i=0; i<L; i++) {//update probability
        for(a=b=j=0; j<N; j++) {
            c = A[i][j>>1] * C[i][j] * B[i][j&(M-1)];
            if(j&1) a+=c;
            else    b+=c;
        }
        if(b+=a) p[i] = a/b;
        else p[i] = 0.5;
    }
}

long ConvolCode::decode(Vec<GF2>& s, const Vec<double>& t,
                        long IMAX, double EPS)
// output: s = decoding of t
// input:
//   t = received code with gaussian noise (prob(s==1))
//   IMAX = max number of iterations for decoding
//   EPS = termination criterion: terminate if
//         min(|prob-0.5|) > 0.5-EPS
// return number of iterations
{
    long i,j,k, L(t.length()/2), N(X.NumCols());
    double T(0.5-EPS);
    Vec<double> p;
    Mat<double> P;

    p.SetLength(L,0.5);
    P.SetDims(L,N);
    for(i=k=0; i<L; i++, k+=2)
        for(j=0; j<N; j++)// initilize transition probability
            P[i][j] = (IsOne(X[0][j]) ? t[k] : 1-t[k])
                     *(IsOne(X[1][j]) ? t[k+1] : 1-t[k+1]);
    for(j=1; j<=IMAX; j++) {
        BCJR(p,P);
        for(i=0; i<p.length(); i++)
            if(fabs(p[i]-0.5) <= T) break;
        if(i==p.length()) break;
    }
    s.SetLength(p.length()-m);
    clear(s);
    for(i=0; i<s.length(); i++)
        if(p[i]>=0.5) set(s[i]);
    return j;
}