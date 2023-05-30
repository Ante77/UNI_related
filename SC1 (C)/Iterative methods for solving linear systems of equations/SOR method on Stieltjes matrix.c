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
integer n=100;
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
    printf("\n %ld . Norma reziduala je: %lf \n", k+1, norma_reziduala);
    uvjet=norma_reziduala/norma;
    printf("\n %ld . Norma reziduala je: %.10lf \n", k+1, uvjet);
    k++;
    

}

printf("\n {");
for(i=0; i<n;i++) printf(" %lf ",x_0[i]);
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

doublereal omega=1.25, tocnost=1e-5;
int counter=0;

int i, k, j;
integer n=100;
FILE *f;
doublereal *A=malloc(n*n*sizeof(doublereal));
 
         
f=fopen("stieltjes_matr.txt","r");
for(i=0;i<n*n;i++)
    fscanf(f,"%lf",A+i);
fclose(f);

char uplo_r='U', uplo_l='L';

doublereal *M=malloc(n*n*sizeof(double));
doublereal *N=malloc(n*n*sizeof(double));

for(i=0; i<n; i++) N[i]=0;
for(i=0; i<n; i++) M[i]=0;

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

char side='L', uplo='L', transa='N', diagg='N';
doublereal alphaa=1.0;

dtrsm_(&side, &uplo, &transa, &diagg, &n, &n, &alphaa, M, &n, N, &n);

sor_konvergencija(N, n);

doublereal *b=malloc(n*sizeof(doublereal));
doublereal *x_T=malloc(n*sizeof(doublereal));
         
for(i=0; i<n; i++){
    x_T[i]=1.0;
}

char trans='N';
doublereal alfa=1.0, beta=0.0;
integer incx=1, incy=1;
dgemv_(&trans, &n, &n, &alfa, A, &n, x_T, &incx, &beta, b, &incy);

doublereal nr_b=dnrm2_(&n, b, &incy);
printf("\n Norma od b: %lf \n", nr_b);

doublereal brojnik=sor_rjesavac (A, b, omega, tocnost);

doublereal nazivnik=dnrm2_(&n, x_T, &incx);

doublereal rjesenje=brojnik/nazivnik;
printf("\n %e \n", rjesenje);

return 0;

}
