#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "clapack.h"
#include "fblaswr.h"
#include <math.h>

integer main(integer argc, char *argv[]){

    integer n=21;
    integer i, j;

    doublereal* A=malloc(n*n*sizeof(doublereal));
    doublereal xi, yi;
    
    for (i=0; i<21; i++){
        for (j=0; j<21; j++){
            xi=-5+0.5*i;
            yi=-5+0.5*j;
            A[i+j*n]=(xi*xi*yi-xi*xi-yi*yi+175)/250;
        }
    }
    
    doublereal* b=malloc(n*sizeof(doublereal));
    b[0]=10;
    for (i=1; i<21; i++){
        b[i]=1;
    }
    
    int* jpvt=malloc(n*sizeof(int));
    doublereal* tau=malloc(n*sizeof(doublereal));
    integer lwork=3*n+1;
    doublereal* work=malloc(lwork*sizeof(doublereal));
    integer info;
    
    dgeqp3_(&n, &n, A, &n, jpvt, tau, work, &lwork, &info);

    doublereal tol=21e-16;
    integer br=0;
    while(1)
         {
            if(abs(A[br+br*n])<tol)  break;    
            br++;
         }
    
    doublereal* R=malloc(br*n*sizeof(doublereal));
    
    for (i=0; i<br; i++){
        for (j=i; j<n; j++){
            R[i+j*br]=A[i+j*n];
        }
    }
    
    integer lwork2=n;
    integer info2;
    
    dorgqr_(&n, &n, &n, A, &n, tau, work, &lwork2, &info2);
    
    char trans='T';
    doublereal alpha=1.0, beta=0.0;
    doublereal* y=malloc(n*sizeof(doublereal));
    integer incx=1;
    
    dgemv_(&trans, &n, &n, &alpha, A, &n, b, &incx, &beta, y, &incx);

    doublereal* b2=malloc(br*sizeof(doublereal));
    
    for (i=0; i<br; i++){
        b2[i]=y[i];
    }
    
    doublereal* tau3=malloc(br*sizeof(doublereal));
    integer lwork3=br*n;
    integer work3=malloc(lwork3*sizeof(doublereal));
    integer info3;
    
    dgelqf_(&br, &n, R, &br, tau3, work3, &lwork3, &info3);
// chill

    char uplo='L', side='L', trans2='N', diag='N';
    integer m=1;
    
    dtrsm_(&side, &uplo, &trans2, &diag, &br, &m, &alpha, R, &br, b2, &br);

    integer lwork4=br*n;
    doublereal* work4=malloc(lwork4*sizeof(doublereal));
    integer info4;
    
    dorglq_(&br, &n, &br, R, &br, tau3, work4, &lwork4, &info4);

    doublereal* v=malloc(n*sizeof(doublereal));
    char trans3='T';
    
    dgemv_(&trans3, &br, &n, &alpha, R, &br, b2, &m, &beta, v, &m);
 
    doublereal* rjesenje=malloc(n*sizeof(doublereal));
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(jpvt[j]==i+1){
                rjesenje[i]=v[j];
            }
        }
    }
    
    printf("\n \n {");
    for(i=0; i<n-1; i++){
        printf(" %lf,", rjesenje[i]);
    }
    printf(" %lf }^T", rjesenje[n-1]);
    printf("\n \n");
 
    return 0;
}      