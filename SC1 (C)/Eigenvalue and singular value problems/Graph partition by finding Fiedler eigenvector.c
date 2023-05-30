#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "clapack.h"
#include "fblaswr.h"
#include <math.h>

integer main(integer argc, char *argv[]){

    integer n=7;

/*
w(1)=9
w(2)=10
w(3)=9
w(4)=14
w(5)=11
w(6)=14
w(7)=9

*/

    integer i, j;

    doublereal L[]={9, -2, -3, -4, 0, 0, 0,
    -2, 10, 0, -7, -1, 0, 0,
    -3, 0, 9, -3, 0, -2, -1,
    -4, -7, -3, 14, 0, 0, 0,
    0, -1, 0, 0, 11, -7, -3,
    0, 0, -2, 0, -7, 14, -5,
    0, 0, -1, 0, -3, -5, 9};

    for(i=0; i<n; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %lf ", L[i+j*n]);
        }
        printf(" \n ");
    }

    char jobz='V', range='I', uplo='U';
    integer il=2, iu=2, m=1;
    doublereal abstol=1e-8;
    doublereal* W=malloc(n*sizeof(doublereal));
    doublereal* Z=malloc(n*sizeof(doublereal));
    integer lwork=8*n;
    doublereal* work=malloc(lwork*sizeof(doublereal));
    integer* iwork=malloc(5*n*sizeof(integer));
    integer* ifail=malloc(n*sizeof(integer));
    integer info;
    doublereal vl, vu;

    dsyevx_(&jobz, &range, &uplo, &n, L, &n, &vl, &vu, &il, &iu, &abstol, &m, W, Z, &n, work, &lwork, iwork, ifail, &info);
    
    printf("\n Svojsveni vektor u je: \n");
    for(i=0; i<n; i++) printf(" %lf ", Z[i]);
    printf("\n \n");
    printf(" za svojstvenu vrijednost %lf ", W[0]);

    doublereal P[]={9, -2, -3, -4, 0, 0, 0,
    -2, 10, 0, -7, -1, 0, 0,
    -3, 0, 9, -3, 0, -2, -1,
    -4, -7, -3, 14, 0, 0, 0,
    0, -1, 0, 0, 11, -7, -3,
    0, 0, -2, 0, -7, 14, -5,
    0, 0, -1, 0, -3, -5, 9};

    doublereal D[]={1/sqrt(9), 0, 0, 0, 0, 0, 0, 0, 1/sqrt(10), 0, 0, 0, 0, 0, 0, 0, 1/sqrt(9), 0, 0, 0, 0, 0, 0, 0, 1/sqrt(14), 0, 0, 0, 0, 0, 0, 0, 1/sqrt(11), 0, 0, 0, 0, 0, 0, 0, 1/sqrt(14), 0, 0, 0, 0, 0, 0, 0, 1/sqrt(9)};

    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        P[i+n*j]=D[i+n*i]*P[i+n*j];
      }
    }
    
    for(j=0;j<n;j++){
      for(i=0;i<n;i++){
        P[i+n*j]=D[j+n*j]*P[i+n*j];
      }
    }

    for(i=0; i<n; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %lf ", P[i+j*n]);
        }
        printf(" \n ");
    }
/*
    doublereal D2[]={1/sqrt(9), 1/sqrt(10), 1/sqrt(9), 1/sqrt(14), 1/sqrt(11), 1/sqrt(14), 1/sqrt(9)};

    char trans='N';
    doublereal alpha=1.0, beta=0.0;
    doublereal* pom2=malloc(n*n*sizeof(doublereal));
    doublereal* pom3=malloc(n*n*sizeof(doublereal));

    dgemm_(&trans, &trans, &n, &n, &n, &alpha, D, &n, P, &n, &beta, pom2, &n);
    dgemm_(&trans, &trans, &n, &n, &n, &alpha, pom2, &n, D2, &n, &beta, pom3, &n);
*/
    doublereal vl2, vu2;
    doublereal* W2=malloc(n*sizeof(doublereal));
    doublereal* Z2=malloc(n*sizeof(doublereal));
    doublereal* work2=malloc(lwork*sizeof(doublereal));
    integer* iwork2=malloc(5*n*sizeof(integer));
    integer* ifail2=malloc(n*sizeof(integer));

    dsyevx_(&jobz, &range, &uplo, &n, P, &n, &vl2, &vu2, &il, &iu, &abstol, &m, W2, Z2, &n, work2, &lwork, iwork2, ifail2, &info);

    printf("\n Svojsveni vektor u je: \n");
    for(i=0; i<n; i++) printf(" %lf ", Z2[i]);
    printf("\n \n");
    printf(" za svojstvenu vrijednost %lf ", W2[0]);
    printf("\n \n");

}