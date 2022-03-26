#include "util.h"

void indexx(const int arr[], int indx[], int n)
// argsort
// W.Press et al., "Numrical Recipes", section 8.4
{
    const int M=7,NSTACK=50;
    int i,indxt,ir,j,k,jstack=-1,l=0;
    int a;
    int istack[NSTACK];

    ir=n-1;
    for (j=0;j<n;j++) indx[j]=j;
    for (;;) {
        if (ir-l < M) {
            for (j=l+1;j<=ir;j++) {
                indxt=indx[j];
                a=arr[indxt];
                for (i=j-1;i>=l;i--) {
                    if (arr[indx[i]] <= a) break;
                    indx[i+1]=indx[i];
                }
                indx[i+1]=indxt;
            }
            if (jstack < 0) break;
            ir=istack[jstack--];
            l=istack[jstack--];
        } else {
            k=(l+ir) >> 1;
            SWAP(indx[k],indx[l+1]);
            if (arr[indx[l]] > arr[indx[ir]]) {
                SWAP(indx[l],indx[ir]);
            }
            if (arr[indx[l+1]] > arr[indx[ir]]) {
                SWAP(indx[l+1],indx[ir]);
            }
            if (arr[indx[l]] > arr[indx[l+1]]) {
                SWAP(indx[l],indx[l+1]);
            }
            i=l+1;
            j=ir;
            indxt=indx[l+1];
            a=arr[indxt];
            for (;;) {
                do i++; while (arr[indx[i]] < a);
                do j--; while (arr[indx[j]] > a);
                if (j < i) break;
                SWAP(indx[i],indx[j]);
            }
            indx[l+1]=indx[j];
            indx[j]=indxt;
            jstack += 2;
            if (jstack >= NSTACK) error("NSTACK too small in indexx.");
            if (ir-i+1 >= j-l) {
                istack[jstack]=ir;
                istack[jstack-1]=i;
                ir=j-1;
            } else {
                istack[jstack]=j-1;
                istack[jstack-1]=l;
                l=i;
            }
        }
    }
}

#include<NTL/ZZ.h>

void indexr(const int arr[], int indx[], int n)
// argsort but random choice from same values
{
    int i,j,k(0);
    indexx(arr, indx, n);
    for(i=1; i<=n; i++) {
        if(i<n && arr[indx[i]] == arr[indx[i-1]]) continue;
        for(j=k; j<i; j++)
            SWAP(indx[j], indx[NTL::RandomBnd(i-j)+j]);
        k=i;
    }
}
