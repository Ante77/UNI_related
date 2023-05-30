#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

/* Za isprobati ostale norme promijeniti character u linijama 13 i 86. Ostali zeleni komentari su nebitni.  */

doublereal sor_norma (doublereal* matrica){

char norm='F';
integer dim=4;
doublereal *work=malloc(dim*dim*sizeof(double));
doublereal rez;

rez=dlange_(&norm, &dim, &dim, matrica, &dim, work);

return rez;

}

integer sor_konvergencija (doublereal* matrica){

doublereal normie;
normie=sor_norma(matrica);
printf("\n Norma matrice je: %lf \n",normie);

if(normie<1) return 1;

return 0;

}

void sor_rjesavac (doublereal* T, doublereal* A, doublereal* b){

char jednou='U';
doublereal omega;
omega=1.05;

integer i, j;
integer n=4;
doublereal pom;

doublereal *x_0=malloc(n*n*sizeof(double));
doublereal *x_1=malloc(n*n*sizeof(double));

for(i=0; i<n; i++) x_0[i]=0;
for(i=0; i<n; i++) x_1[i]=0;

dlacpy_(&jednou, &n, &n, x_0, &n, x_1, &n);

for(i=0; i<n; i++){
        
    x_1[i]*=(1-omega);
    pom=b[i];
        
    for(j=0; j<i; j++){
        pom-=A[i+j*n]*x_1[j];
    }
    for(j=i+1; j<n; j++){
        pom-=A[i+j*n]*x_1[j];
    }
        
    x_1[i]+=(pom*omega)/A[i+i*n];

}

doublereal kriterij;
doublereal vekt_norma=0, mat_norma=0;

char norm='F';
integer dim=1;
integer dim_matr=4;
doublereal *work=malloc(dim_matr*dim_matr*sizeof(double));

vekt_norma=dlange_(&norm, &dim_matr, &dim, x_1, &dim_matr, work);
mat_norma=sor_norma(T);

kriterij=(log((0.00001*(1-mat_norma))/(vekt_norma)))/(log(mat_norma));

integer k=1;

while(k<kriterij){

    k++;
    for(i=0; i<n; i++){
        
        x_0[i]*=(1-omega);
        pom=b[i];
        
        for(j=0; j<i; j++){
            pom-=A[i+j*n]*x_0[j];
        }
        for(j=i+1; j<n; j++){
            pom-=A[i+j*n]*x_0[j];
        }
        
        x_0[i]+=(pom*omega)/A[i+i*n];

    }

}

printf("Matricna norma je: %lf \n \n", mat_norma);

printf("Rjesenje je: ");
printf("{");
for(i=0; i<n; i++)
{
printf(" %lf ", x_0[i]);
}
printf("}\n \n");
printf("Broj iteracija: %ld\n \n", k);

}

main(integer argc, char*argv[]){

integer n=4;
integer dimenzija=16;

doublereal *A=malloc(n*n*sizeof(double));
doublereal *B=malloc(n*sizeof(double));
A[0]=101.0; A[1]=-4.0; A[2]=8.0; A[3]=12.0; A[4]=-4.0; A[5]=20.0; A[6]=-7.0; A[7]=3.0;
A[8]=8.0; A[9]=-7.0; A[10]=78.0; A[11]=32.0; A[12]=12.0; A[13]=3.0; A[14]=32.0; A[15]=113.0;   
B[0]=117.0; B[1]=12.0; B[2]=111.0; B[3]=160.0;   

doublereal omega=1.05;
doublereal tocnost=0.00001;

char uplo_r='U', uplo_l='L';

doublereal *D=malloc(n*n*sizeof(double));
doublereal *M=malloc(n*n*sizeof(double));
doublereal *N=malloc(n*n*sizeof(double));
integer i, j;

for(i=0; i<n; i++) N[i]=0;
for(i=0; i<n; i++) D[i]=0;
for(i=0; i<n; i++) M[i]=0;

dlacpy_(&uplo_r, &n, &n, A, &n, N, &n);

dlacpy_(&uplo_l, &n, &n, A, &n, M, &n);

integer k;

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

printf("\n Matrica T je: \n");

dtrsm_(&side, &uplo, &transa, &diagg, &n, &n, &alphaa, M, &n, N, &n);

for(i=0; i<n; i++){
    for(j=0; j<n; j++){
        printf(" %lf ",N[i+j*n]);
    }
    printf("\n");
}

printf("\n Vektor B je: \n"); 

for(i=0; i<n; ++i)
{
    printf(" %lf ",B[i]);
}
printf("\n \n");

sor_rjesavac(N, A, B);

return 0;

}