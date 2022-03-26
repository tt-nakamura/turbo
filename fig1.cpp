#include "ConvolCode.h"
#include "TurboCode.h"
#include "LDPCCode.h"
#include<fstream>

main() {
    long N(512),K(256),IMAX(8);
    long n(8),M(256);
    double s1(0.6), s2(1);
    double sigma, e[3], w, ds((s2-s1)/n);
    long i,j,k;
    Vec<GF2> u,v;
    Vec<double> y;
    std::ofstream f("fig1.txt");

    ConvolCode c;
    TurboCode b(K);
    LDPCCode l(N,K);

    for(i=0; i<=n; i++) {
        sigma = s1 + i*ds;
        std::cout << sigma;
        f << sigma;
        for(k=0; k<3; k++) e[k] = 0;
        for(j=0; j<M; j++) {
            random(u,K);
            // Convol
            c.encode(v,u);
            AddNoise(y,v,sigma);
            c.decode(v,y,IMAX);
            w = weight(v-=u);
            e[0] += w/K;
            // Turbo
            b.encode(v,u);
            AddNoise(y,v,sigma);
            b.decode(v,y,IMAX);
            w = weight(v-=u);
            e[1] += w/K;
            // LDPC
            l.encode(v,u);
            AddNoise(y,v,sigma);
            l.decode(v,y,IMAX);
            w = weight(v-=u);
            e[2] += w/K;
        }
        for(k=0; k<3; k++) {
            e[k] /= M;
            std::cout << ' ' << e[k];
            f << ' ' << e[k];
        }
        std::cout << std::endl;
        f << std::endl;
    }
}