#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "clapack.h"
#include "fblaswr.h"
#include <math.h>

integer main(integer argc, char *argv[]){

    integer m=10, n=2;
    integer i, j;

    doublereal A[]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    doublereal B[]={3.5, 4.9, 6.8, 9.3, 10.9, 13.4, 15.1, 16.7, 19.0, 21.2};

    char jobu='A', jobvt='A';
    doublereal* S=malloc(n*sizeof(doublereal));
    doublereal* U=malloc(m*m*sizeof(doublereal));
    doublereal* Vt=malloc(n*n*sizeof(doublereal));
    integer lwork=3*n+m;
    doublereal* work=malloc(lwork*sizeof(doublereal));
    integer info; 

    dgesvd_(&jobu, &jobvt, &m, &n, A, &m, S, U, &m, Vt, &n, work, &lwork, &info);

    printf("Prva dva stupca od U:\n");
    for(i=0; i<m; i++){
        printf("\n");
        for(j=0; j<2; j++){
            printf(" %lf ", U[i+j*m]);
        }
    }
    printf("\n \n");
    printf("Dijagonala od S:\n \n");
    for(i=0; i<n; i++){
        printf(" %lf ", S[i]);
    }
    printf("\n");

    printf("\n \n");
    
    doublereal* V=malloc(n*n*sizeof(doublereal));

    printf("Matrica V je:\n \n");
    for(i=0; i<n; i++){
        printf("\n");
        for(j=0; j<n; j++){
           // printf(" %lf ", Vt[j+i*n]);
            V[i+j*n]=Vt[j+i*n];
        }
    }
    printf("\n");

    printf("Matrica V je:\n \n");
    for(i=0; i<n; i++){
        printf("\n");
        for(j=0; j<n; j++){
            printf(" %lf ", V[i+j*n]);
        }
    }
    printf("\n");

    doublereal* sigma=malloc(n*n*sizeof(doublereal));
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            sigma[i+j*n]=0.0;
        }
        sigma[i+i*n]=1/S[i];
    }

    doublereal* Upd=malloc(m*2*sizeof(doublereal));
    for(i=0; i<m*2; i++){
        Upd[i]=U[i];
    }

    doublereal alpha=1.0, beta=0.0;

    char transa='T', transb='N';
    integer k=1;
    doublereal* rez1=malloc(n*k*sizeof(doublereal));

// prvo mnozenje
    dgemm_(&transa, &transb, &n, &k, &m, &alpha, Upd, &m, B, &m, &beta, rez1, &n);

    printf("Prvo mnozenje je:\n \n");
    for(i=0; i<n; i++){
        printf("\n");
        for(j=0; j<k; j++){
            printf(" %lf ", rez1[i+j*n]);
        }
    }
    printf("\n");
// drugo mnozenje
    doublereal* rez2=malloc(n*k*sizeof(doublereal));
    dgemm_(&transb, &transb, &n, &k, &n, &alpha, sigma, &n, rez1, &n, &beta, rez2, &n);

    printf("Drugo mnozenje je:\n \n");
    for(i=0; i<n; i++){
        printf("\n");
        for(j=0; j<k; j++){
            printf(" %lf ", rez2[i+j*n]);
        }
    }
    printf("\n");
// trece mnozenje
    doublereal* rez3=malloc(n*k*sizeof(doublereal));
    dgemm_(&transb, &transb, &n, &k, &n, &alpha, V, &n, rez2, &n, &beta, rez3, &n);

    printf("\n \n");
    printf("Rezultat:\n \n");
    for(i=0; i<n; i++){
        printf(" %lf ", rez3[i]);
    }
    printf("\n");

}