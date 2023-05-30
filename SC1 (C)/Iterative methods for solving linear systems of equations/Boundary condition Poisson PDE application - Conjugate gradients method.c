#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

// qr funkcija izracuna gornjetr R ali ne i Q nego sprema informacije o HH reflektorima pa dormqr koristimo da dobijemo matricu A
// u dormqr tau i a su iz qr funkcije a c je nasa dijagonalna matrica

void cg(integer n, doublereal* A, doublereal* b, doublereal* x_0, doublereal tol){
    
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
// cg algoritam    
doublereal *d=malloc(n*sizeof(doublereal));
doublereal *adk=malloc(n*sizeof(doublereal));
doublereal *rezerva=malloc(n*sizeof(doublereal));
integer k=0;
char uplo='U';
doublereal alfica=1.0, betica=0.0;
integer incx=1, incy=1;

for(i=0; i<n; i++){
    d[i]=0.0;
    adk[i]=0.0;
    rezerva[i]=0.0;
}

dcopy_(&n, rezidual, &incb, d, &incr);

doublereal stari, produkt_dkaiAdka, novi;
doublereal aleph_k, minus_aleph_k, betath_k;

stari=ddot_(&n, rezidual, &incx, rezidual, &incy);

while(uvjet>tol){
/*
    printf("\n x_0 je: \n");
    for(i=0; i<n; i++){
        printf(" %lf ", x_0[i]);
    }
*/        
    dsymv_(&uplo, &n, &alfica, A, &n, d, &incx, &betica, adk, &incy); // produkt A*d_k spremljen u adk
    
    produkt_dkaiAdka=ddot_(&n, d, &incx, adk, &incy);
    
//    printf("\n Produkt r_k ova je %e \n", produkt_rkova);
//    printf("\n Produkt d_k i Ad_k je %e \n", produkt_dkaiAdka);
    
    aleph_k=stari/produkt_dkaiAdka;
//    printf("\n Aleph k je %e \n", aleph_k);
    minus_aleph_k=aleph_k*(-1);
//    printf("\n Minus aleph k je %e \n", minus_aleph_k);
    
    daxpy_(&n, &aleph_k, d, &incx, x_0, &incy);  // racunanje x_(k+1)
    daxpy_(&n, &minus_aleph_k, adk, &incx, rezidual, &incy);  // racunanje r_(k+1)
    
    novi=ddot_(&n, rezidual, &incx, rezidual, &incy);
    betath_k=novi/stari;
    
//    printf("\n Betica k je %e \n", betath_k);

    dcopy_(&n, rezidual, &incx, rezerva, &incy);
    daxpy_(&n, &betath_k, d, &incx, rezerva, &incy); // pravim d_(k+1)
    dcopy_(&n, rezerva, &incx, d, &incy);
    
    stari=novi;
    norma_reziduala=sqrt(novi);
    uvjet=norma_reziduala/norma_b;
    
    printf("\n Relativna norma reziduala u %ld koraku je: %e \n", k+1, uvjet);
    
    k++;
}

printf("\n Rjesenje je: \n");
printf("\n {");
for(i=0; i<9; i++){
  for(j=0; j<9; j++){
    printf(" %lf ", x_0[i+j*9]);
  }
  printf("\n");
}
printf("} \n");

printf("Broj iteracija: %d\n \n", k);

}

int main(integer argc, char*argv[]){

integer n=81, i, j;

doublereal tol=1e-8;

// scanf("tol = %lf", &tol);

doublereal *A=malloc(n*n*sizeof(double));
doublereal *B=malloc(n*sizeof(double));
doublereal *x_0=malloc(n*sizeof(double));

for(i=0; i<n; i++) x_0[i]=0.0;

// Pravim matricu A

for(i=0; i<n; i++){
  for(j=0; j<n; j++){
    if(i==j){
      A[i+j*n]=4*100;
      if(i!=0 && i%9!=0) A[i+j*n-1]=(-1)*100;
      if(i!=n-1 && i%9!=8) A[i+j*n+1]=(-1)*100;
    }
    if(i==j+9) A[i+j*n]=(-1)*100;
    if(j==i+9) A[i+j*n]=(-1)*100;
  }
}
for(i=0; i<n; i++){
  B[i]=0;
  B[40]=10000;
}

cg(n, A, B, x_0, tol);

return 0;
}