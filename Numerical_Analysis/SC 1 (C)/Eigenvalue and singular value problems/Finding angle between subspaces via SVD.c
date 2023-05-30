#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "clapack.h"
#include "fblaswr.h"
#include <math.h>

integer main(integer argc, char *argv[]){

    integer m=4, n=2;
    integer i, j;

    doublereal A[]={1, 0, 0, 0, 1, 1, 1, 1};
    doublereal B[]={1, -1, 1, -1, 0, 1, 0, 1};

    // TREBAMO QR FAKTORIZACIJU A I B MATRICA

    doublereal* tauA=malloc(n*sizeof(doublereal));
    doublereal* tauB=malloc(n*sizeof(doublereal));
    integer lwork5=n*n;
    integer lwork6=n*n;
    doublereal* workA=malloc(lwork5*sizeof(doublereal));
    doublereal* workB=malloc(lwork6*sizeof(doublereal));
    integer info1;
    integer info2;
    
    dgeqrf_(&m, &n, A, &m, tauA, workA, &lwork5, &info1);
    dgeqrf_(&m, &n, B, &m, tauB, workB, &lwork6, &info2);

    doublereal* workB2=malloc(lwork5*sizeof(doublereal));
    integer info3;
    doublereal* workA2=malloc(lwork6*sizeof(doublereal));
    integer info4;


    dorgqr_(&m, &n, &n, B, &m, tauB, workB2, &lwork5, &info3);
    dorgqr_(&m, &n, &n, A, &m, tauA, workA2, &lwork6, &info4);

    printf("Matrica Qb je: \n");
    for(i=0; i<8; i++){
            printf(" %lf ", B[i]);
    }
    printf("Matrica Qa je: \n");
    for(i=0; i<8; i++){
            printf(" %lf ", A[i]);
    }


   /* char side='L', trans='T';
    integer lwork3=m;
    doublereal* workA2=malloc(lwork3*sizeof(doublereal));
    integer info4;
*/
   // dormqr_(&side, &trans, &n, &m, &n, A, &m, tauA, B, &n, workA2, &lwork3, &info4);

    // C je 0.5 0.5
    //      -0.2887 0.8660

    doublereal* C=malloc(n*n*sizeof(doublereal));

    char transa='T', transb='N';
    doublereal alpha=1.0, beta=0.0;

    dgemm_(&transa, &transb, &n, &n, &m, &alpha, A, &m, B, &m, &beta, C, &n);

   // C[0]=0.5; C[2]=0.5;

    printf("Matrica C je: \n");
    for(i=0; i<n; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %lf ", C[i+j*n]);
        }
        printf(" \n ");
    }

    // tu smo, B=Q_a^t*Q_b

    char jobu='A', jobvt='A';

    doublereal* S=malloc(n*sizeof(doublereal));
    doublereal* U=malloc(n*n*sizeof(doublereal));
    doublereal* Vt=malloc(n*n*sizeof(doublereal));
    integer lwork2=5*n;
    doublereal* work2=malloc(lwork2*sizeof(doublereal));
    integer info; 

    dgesvd_(&jobu, &jobvt, &n, &n, C, &n, S, U, &n, Vt, &n, work2, &lwork2, &info);

    printf("\n \n");
    for(i=0; i<n; i++){
        printf(" %e ", S[i]);
    }
    printf("\n \n");

  //  printf("\n Infosi: %ld %ld %ld %ld %ld \n", info1, info2, info3, info4, info);

    doublereal* X=malloc(m*n*sizeof(doublereal));
    doublereal* Y=malloc(m*n*sizeof(doublereal));

    dgemm_(&transb, &transb, &m, &n, &n, &alpha, A, &m, U, &n, &beta, X, &m);
    dgemm_(&transb, &transa, &m, &n, &n, &alpha, B, &m, Vt, &n, &beta, Y, &m);

   /* printf("Matrica X je: \n");
    for(i=0; i<m; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %lf ", X[i+j*n]);
        }
        printf(" \n ");
    }

    printf("Matrica Y je: \n");
    for(i=0; i<m; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %lf ", Y[i+j*n]);
        }
        printf(" \n ");
    }
    */

    printf("\n \n");
    for(i=0; i<m*n; i++){
        printf(" %lf ", X[i]);
    }

    printf("\n \n");
    for(i=0; i<m*n; i++){
        printf(" %lf ", Y[i]);
    }
    
    printf("\n \n Kutevi su: %lf %lf \n", acos(S[0]), acos(S[1]));

}