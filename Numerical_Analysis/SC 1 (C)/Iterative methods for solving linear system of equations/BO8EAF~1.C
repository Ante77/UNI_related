#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

doublereal sor_norma (doublereal* matrica, integer n){

char norm='F';
doublereal *work=malloc(n*n*sizeof(double));
doublereal rez;

rez=dlange_(&norm, &n, &n, matrica, &n, work);

printf("\n Norma matrice je: %lf \n",rez);

return rez;

}

integer sor_konvergencija (doublereal* matrica, integer n){

doublereal normie=sor_norma(matrica, 100);
printf("\n Nova norma matrice je: %lf \n",normie);

if(normie<1) return 1;

return 0;

}

// RJESAVAC 

doublereal sor_rjesavac (doublereal* A, doublereal* b, doublereal omega, doublereal tocnost){

integer i, j;
integer n=99;
doublereal pom;

doublereal *x_0=malloc(n*sizeof(doublereal));
doublereal *r=malloc(n*sizeof(doublereal));
doublereal *y=malloc(n*sizeof(doublereal));
doublereal *ikste=malloc(n*sizeof(doublereal));
         
for(i=0; i<n; i++){
    x_0[i]=0.0;
    y[i]=0.0;
    ikste[i]=1.0;
}

integer incb=1, incr=1;

dcopy_(&n, b, &incb, r, &incr);

char trans='n';
doublereal alfa=-1.0, beta= 1.0;
dgemv_(&trans, &n, &n, &alfa, A, &n, x_0, &incb, &beta, b, &incb);

doublereal norma=dnrm2_(&n, r, &incr);
doublereal norma_reziduala=dnrm2_(&n, b, &incb);
doublereal uvjet=norma_reziduala/norma;

printf("\n Norma reziduala je: %lf \n", norma_reziduala);
printf("\n Norma b je: %lf \n", norma);

printf("\n Uvjet je: %lf \n", uvjet);

integer k=0;

while(uvjet>tocnost){

    for(i=0; i<n; i++){
        
        x_0[i]*=(1-omega);
        pom=r[i];
        
        for(j=0; j<i; j++){
            pom-=A[i+j*n]*x_0[j];
        }
        for(j=i+1; j<n; j++){
            pom-=A[i+j*n]*x_0[j];
        }
        
        x_0[i]+=(pom*omega)/A[i+i*n];

    }
    
    dcopy_(&n, r, &incb, y, &incr);
    dgemv_(&trans, &n, &n, &alfa, A, &n, x_0, &incb, &beta, y, &incb);
    norma_reziduala=dnrm2_(&n, y, &incb);
    printf("\n %ld . Norma reziduala je: %e \n", k+1, norma_reziduala);
    uvjet=norma_reziduala/norma;
    printf("\n %ld . Uvjet  je: %e \n", k+1, uvjet);
    k++;
    

}

doublereal *vektor_razlike=malloc(n*sizeof(doublereal));
for(i=1; i<n+1; i++){
    vektor_razlike[i-1]=0.01*i*cos(0.01*i)-x_0[i-1];
}

printf("\n Razlika aproksimacije i rjesenja je: \n");
printf("\n {");
for(i=0; i<n;i++) printf(" %e ",vektor_razlike[i]);
printf("} \n");

printf("Broj iteracija: %d\n \n", k);

doublereal *p=malloc(n*sizeof(doublereal));
         
for(i=0; i<n; i++){
    p[i]=ikste[i]-x_0[i];
}

doublereal brojnick=dnrm2_(&n,p,&incb);

return brojnick;

}

// NJU MAIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

int main(integer argc, char*argv[]){

doublereal omega=1.0, tocnost=1e-8;
int counter=0;

int i, k, j;
integer n=99;

doublereal *A=malloc(n*n*sizeof(double));
doublereal *B=malloc(n*sizeof(double));
for(i=0; i<n*n; i++){
    A[i]=0.0;
}
for(i=0; i<n; i++){
    A[i+i*n]=1.9999;
}
for(i=0; i<n-1; i++){
    A[i+i*n+1]=-1;
}
for(i=1; i<n; i++){
    A[i+i*n-1]=-1;
}

for(i=1; i<n; i++){
    B[i-1]=0.0002*sin(0.01*i);
}
B[n-1]=0.0002*sin(0.99)+cos(1);

integer infosboi, lda=n;
char uplo='U';
/*
dpotrf_(&uplo, &n, A, &lda, &infosboi);
if(infosboi) printf("\n MATRICA NIJE SIMETRICNA OMG!!! \n");
if(!infosboi) printf("\n Matrica je simetricna ako je za vjerovat ovim robotima. \n");
*/
printf("\n Matrica A je: \n");
for(i=0; i<n; i++){
    for(j=0; j<n; j++){
        printf(" %lf ",A[i+j*n]);
    }
    printf("\n");
}
printf("\n Vektor B je: \n");
for(i=0; i<n; i++){
    printf(" %e ",B[i]);
    
}
    

doublereal *M=malloc(n*n*sizeof(double));
doublereal *N=malloc(n*n*sizeof(double));

for(i=0; i<n; i++) N[i]=0;
for(i=0; i<n; i++) M[i]=0;

char uplo_r='U', uplo_l='L';

dlacpy_(&uplo_r, &n, &n, A, &n, N, &n);

dlacpy_(&uplo_l, &n, &n, A, &n, M, &n);

for(i=0; i<n; i++){
    for(j=0; j<n; j++){
        if(i<j) M[i*n+j]*=omega;
    }
}

for(i=0; i<n; i++){
    for(j=0; j<n; j++){
        if(i==j){
            N[i+j*n]*=(1-omega);
        }
        if(i<j){
            N[i+j*n]=(-omega)*N[i+j*n];
        }
    }
}

char side='L', uplol='L', transa='N', diagg='N';
doublereal alphaa=1.0;

dtrsm_(&side, &uplol, &transa, &diagg, &n, &n, &alphaa, M, &n, N, &n);

sor_konvergencija(N, n);

/*
char trans='N';
doublereal alfa=1.0, beta=0.0;
integer incx=1, incy=1;
dgemv_(&trans, &n, &n, &alfa, A, &n, x_T, &incx, &beta, b, &incy);
*/
integer incx=1, incy=1;
doublereal nr_b=dnrm2_(&n, B, &incy);
printf("\n Norma od b: %lf \n", nr_b);

doublereal brojnik=sor_rjesavac (A, B, omega, tocnost);
/*
doublereal nazivnik=dnrm2_(&n, x_T, &incx);

doublereal rjesenje=brojnik/nazivnik;
*/
return 0;

}