#include <stdio.h>
#include <stdlib.h>
#include "f2c.h"
#include "clapack.h"
#include "fblaswr.h"
#include <math.h>

doublereal SodA (doublereal* A, integer n){
    doublereal sum=0, rez;
    integer i, j;
    
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i!=j){
                sum+=A[i+j*n]*A[i+j*n];
            }
        }
    }
    rez=sqrt(sum);
    return rez;
}

void jacobisd (integer n, doublereal* A, doublereal tol){

    char norm='F';
    doublereal* work=malloc(n*n*sizeof(doublereal));
    
    doublereal norma=dlange_(&norm, &n, &n, A, &n, work);
    doublereal uvjet=norma*tol;

    integer p, q, k, i=0, br1, br2;

    doublereal S=SodA(A,n);

    printf("\n S(A)= %e u koraku %ld \n", S, i);
    
    doublereal tau, t, c, s, app, aqq, apq, pom;
    integer sgn;

    while(S>uvjet){

        for(p=0; p<n-1; p++){
            for(q=p+1; q<n; q++){
                if(A[p+q*n]!=0){
                    tau=(A[q+q*n]-A[p+p*n])/(2*A[p+q*n]);
                    
                    if(tau>=0) sgn=1;
                    if(tau<0) sgn=-1;
                    t=sgn/(fabs(tau)+sqrt(1+pow(tau,2)));
                    c=1/(sqrt(1+pow(t,2)));
                    s=t*c;
                    }
                    app=A[p+p*n]; apq=A[p+q*n]; aqq=A[q+q*n];
                    app=app-t*apq;
                    aqq=aqq+t*apq;
                    for(k=0; k<n; k++){
                        pom=A[k+p*n];
                        A[k+p*n]=c*pom-s*A[k+q*n];
                        A[k+q*n]=s*pom+c*A[k+q*n];
                        A[p+k*n]=A[k+p*n];
                        A[q+k*n]=A[k+q*n];
                    }
                    A[p+q*n]=0; A[q+p*n]=0;
                    A[p+p*n]=app; A[q+q*n]=aqq;
            }
        }

        S=SodA(A,n);
        i++;
        printf("\n S(A)= %e u koraku %ld \n", S, i);

    }
    
    for(br1=0; br1<n; br1++){
        printf(" \n ");
        for(br2=0; br2<n; br2++){
            printf(" %e ", A[br1+br2*n]);
        }
        printf(" \n ");
    }

}

integer main(integer argc, char *argv[]){

    doublereal tol=4e-16;
   // scanf("Unesi tol: %lf", &tol);

    integer n=4;
    integer i, j;

    doublereal* U=malloc(n*n*sizeof(doublereal));
    doublereal D[]={-10, 0, 0, 0, 0, -5, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0.2}; 
    
    integer idist=3, iseed[4]={125, 3698 , 98 , 11};
    integer nxn=n*n;
    dlarnv_(&idist, iseed, &nxn, U);
    
    doublereal* tau=malloc(n*n*sizeof(doublereal));
    doublereal* work=malloc(n*n*sizeof(doublereal));
    integer lwork=n, info;
    
    dgeqrf_(&n, &n, U, &n, tau, work, &lwork, &info);
    dorgqr_(&n, &n, &n, U, &n, tau, work, &n, &info);

    doublereal* pom=malloc(n*n*sizeof(doublereal));
    doublereal* A=malloc(n*n*sizeof(doublereal));

    doublereal alpha=1.0, beta=0.0;
    char trans1='N',trans2='T';
    
    dgemm_(&trans1, &trans1, &n, &n, &n, &alpha, U, &n, D, &n, &beta, pom, &n);
    dgemm_(&trans1, &trans2, &n, &n, &n, &alpha, pom, &n, U, &n, &beta, A, &n);

    for(i=0; i<n; i++){
        printf(" \n ");
        for(j=0; j<n; j++){
            printf(" %e ", A[i+j*n]);
        }
        printf(" \n ");
    }

    jacobisd(n, A, tol);

}