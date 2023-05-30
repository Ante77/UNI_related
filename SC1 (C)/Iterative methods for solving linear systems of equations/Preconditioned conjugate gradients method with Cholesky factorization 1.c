#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

void ic (integer n, doublereal *a){

    int i,j,k;
    for (i=0; i<n; i++){
        for (k=0; k<i; k++)
            a[i+i*n]-=a[k+i*n]*a[k+i*n];
        a[i+i*n]=sqrt(a[i+i*n]);
        for (j=i+1; j<n; j++){
            if (a[i+j*n]!=0){
                for (k=0; k<i; k++)
                    a[i+j*n]-=a[k+i*n]*a[k+j*n];
                a[i+j*n]/=a[i+i*n];
            }
            
        }
    }

}

void pcg(integer n, doublereal* A, doublereal* b, doublereal* x_0, doublereal tol){
    
doublereal *R=malloc(n*n*sizeof(doublereal));
integer incx=1, incy=1;
integer dim=n*n;

dcopy_(&dim, A, &incx, R, &incy);

ic(n,R);

integer i, j;
// pravim rezidual kao u proslom zadatku
integer incb=1, incr=1;
char trans='N';
doublereal alfa=-1.0, beta= 1.0;
doublereal *rezidual=malloc(n*sizeof(doublereal));
dcopy_(&n, b, &incb, rezidual, &incr);
    
dgemv_(&trans, &n, &n, &alfa, A, &n, x_0, &incb, &beta, rezidual, &incr);
    
doublereal norma_reziduala=dnrm2_(&n, rezidual, &incr);
doublereal norma_b=dnrm2_(&n, b, &incb);
doublereal uvjet=norma_reziduala/norma_b;

printf("\n %lf %lf \n", norma_reziduala, norma_b);
// prvi put Mp_0=r_o
doublereal *p=malloc(n*n*sizeof(doublereal));

char uplo1='U', trans1='T', diag1='N';
char uplo2='U', trans2='N', diag2='N';

dcopy_(&n, rezidual, &incb, p, &incr);

dtrsv_(&uplo1, &trans1, &diag1, &n ,R, &n, p, &incx);
dtrsv_(&uplo2, &trans2, &diag2, &n ,R, &n, p, &incx);

// pcg algoritam    
doublereal *d=malloc(n*sizeof(doublereal));
doublereal *adk=malloc(n*sizeof(doublereal));
doublereal *rezerva=malloc(n*sizeof(doublereal));
integer k=0;
char uplo='U';
doublereal alfica=1.0, betica=0.0;

for(i=0; i<n; i++){
    d[i]=0.0;
    adk[i]=0.0;
    rezerva[i]=0.0;
}

dcopy_(&n, p, &incb, d, &incr);

doublereal stari, produkt_dkaiAdka, novi;
doublereal aleph_k, minus_aleph_k, betath_k;

stari=ddot_(&n, rezidual, &incx, p, &incy);

while(uvjet>tol){
/*
    printf("\n x_0 je: \n");
    for(i=0; i<n; i++){
        printf(" %lf ", x_0[i]);
    }
*/        
    dsymv_(&uplo, &n, &alfica, A, &n, d, &incx, &betica, adk, &incy); // produkt A*d_k spremljen u adk
    
    produkt_dkaiAdka=ddot_(&n, d, &incx, adk, &incy);
    
//    printf("\n Stari je %e \n", stari);
//    printf("\n Produkt d_k i Ad_k je %e \n", produkt_dkaiAdka);
    
    aleph_k=stari/produkt_dkaiAdka;
//    printf("\n Aleph k je %e \n", aleph_k);
    minus_aleph_k=aleph_k*(-1);
//    printf("\n Minus aleph k je %e \n", minus_aleph_k);
    
    daxpy_(&n, &aleph_k, d, &incx, x_0, &incy);  // racunanje x_(k+1)
    daxpy_(&n, &minus_aleph_k, adk, &incx, rezidual, &incy);  // racunanje r_(k+1)
    
    dcopy_(&n, rezidual, &incb, p, &incr);
    
    dtrsv_(&uplo1, &trans1, &diag1, &n ,R, &n, p, &incx);
    dtrsv_(&uplo2, &trans2, &diag2, &n ,R, &n, p, &incx);
    
    novi=ddot_(&n, rezidual, &incx, p, &incy);
    betath_k=novi/stari;
    
//    printf("\n Betica k je %e \n", betath_k);

    dcopy_(&n, p, &incx, rezerva, &incy);
    daxpy_(&n, &betath_k, d, &incx, rezerva, &incy); // pravim d_(k+1)
    dcopy_(&n, rezerva, &incx, d, &incy);
    
    stari=novi;
    norma_reziduala=dnrm2_(&n, rezidual, &incr);
    uvjet=norma_reziduala/norma_b;
    
    printf("\n Relativna norma reziduala u %ld koraku je: %e \n", k+1, uvjet);
    
    k++;
}

printf("\n Broj iteracija: %ld \n", k);

printf("\n Rezultat: \n");
for(i=0; i<n; i++){
    printf(" %lf ", x_0[i]);

}

}

int main(integer argc, char*argv[]){

integer n=100, i;
FILE *f;
doublereal *A=malloc(n*n*sizeof(doublereal));

f=fopen("stieltjes_matr.txt","r");
for(i=0;i<n*n;i++)
    fscanf(f,"%lf",A+i);
fclose(f);

doublereal *x_0=malloc(n*n*sizeof(doublereal));
doublereal *x_T=malloc(n*n*sizeof(doublereal));
doublereal *b=malloc(n*n*sizeof(doublereal));
for(i=0; i<n; i++){
    x_0[i]=0.0;
    x_T[i]=1.0;
    b[i]=0.0;
}

char trans='N';
doublereal alfa=1.0, beta=0.0;
integer incx=1, incy=1;

// doublereal tol=0.00000001;
doublereal tol;
scanf("%lf", &tol);

dgemv_(&trans, &n, &n, &alfa, A, &n, x_T, &incx, &beta, b, &incy);
/*
for(i=0; i<n; i++){
    printf(" %lf ", b[i]);
}
*/
pcg(n, A, b, x_0, tol);

return 0;

}