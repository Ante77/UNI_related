#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "clapack.h"
#include "fblaswr.h"
#include <math.h>

integer main(integer argc, char *argv[]){

    doublereal A[]={1.2, 2.9, 5.2, 6.8, 2.1, 4.3, 6.1, 8.1};
    doublereal B[]={1, 3, 5, 7, 2, 4, 6, 8};

    integer m=4, n=2;
    integer i, j;

    // C=A^t*B
    // SVD U^t*C*V=S
    // Q=U*V^t

    char transa='T', transb='N';
    doublereal alpha=1.0, beta=0.0;
    doublereal* C=malloc(n*n*sizeof(doublereal));
    
    printf("\n Kaj te muci njofra? \n");
    dgemm_(&transa, &transb, &n, &n, &m, &alpha, A, &m, B, &m, &beta, C, &n);

    printf("Matrica C je: \n");
    for(i=0; i<n; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %lf ", C[i+j*n]);
        }
        printf(" \n ");
    }
    doublereal* fejkA=malloc(m*n*sizeof(doublereal));
    for(i=0; i<m*n; i++){
        fejkA[i]=A[i];
    }

    char jobu='A', jobvt='A';

    doublereal* S=malloc(n*n*sizeof(doublereal));
    doublereal* U=malloc(n*n*sizeof(doublereal));
    doublereal* Vt=malloc(n*n*sizeof(doublereal));
    integer lwork=5*n;
    doublereal* work=malloc(lwork*sizeof(doublereal));
    integer ldu=n, info; 

    dgesvd_(&jobu, &jobvt, &n, &n, C, &n, S, U, &n, Vt, &n, work, &lwork, &info);

    doublereal* Q=malloc(2*2*sizeof(doublereal));
    char trans='N';

    printf("\n Kaj te muci njofra? \n");
    dgemm_(&trans, &trans, &n, &n, &n, &alpha, U, &n, Vt, &n, &beta, Q, &n);

    printf("\n Matrica Qje: \n");
    for(i=0; i<n; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %lf ", Q[i+j*n]);
        }
        printf(" \n ");
    }

    // AQ-B racunam

    doublereal alfica=1.0, betica=-1.0;

    dgemm_(&trans, &trans, &m, &n, &n, &alfica, fejkA, &m, Q, &n, &betica, B, &m);

    printf("\n Matrica AQ-B je: \n");
    for(i=0; i<n; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %lf ", B[i+j*n]);
        }
        printf(" \n ");
    }

    char norm='F';
    doublereal* work2=malloc(m*m*sizeof(doublereal));

    doublereal rez=dlange_(&norm, &m, &n, B, &m, work2);
    
    printf("\n Norma je: %lf \n", rez);

}