#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

doublereal sor_rjesavac (doublereal* A, doublereal* b, doublereal* x_0, integer n, doublereal omega, integer br_iter){

integer i, j, iter;
doublereal pom;

doublereal *r=malloc(n*sizeof(doublereal));
doublereal *y=malloc(n*sizeof(doublereal));
         
for(i=0; i<n; i++){
    y[i]=x_0[i];
}

integer incb=1, incr=1;
dcopy_(&n, b, &incb, r, &incr);

for(iter=0; iter<br_iter; iter++){
    for(i=0; i<n; i++){
        
        y[i]*=(1-omega);
        pom=r[i];

        for(j=0; j<i; j++){
            pom-=A[i+j*n]*y[j];
        }
        for(j=i+1; j<n; j++){
            pom-=A[i+j*n]*y[j];
        }
        y[i]+=(pom*omega)/A[i+i*n];
    }   
}

doublereal n_inf, max=0;
for(i=0; i<n; i++){
    n_inf=abs(y[i]-1);
    if(n_inf>max){
        max=n_inf;
    }
}

doublereal jedan_kroz_N=1.0/br_iter;
doublereal p=pow(max,jedan_kroz_N);
return max;
}

int main(integer argc, char*argv[]){

integer n=256;
integer p=sqrt(n);
integer br_iter=100;
integer br_iter2=100;
integer i, j;
doublereal* A=malloc(n*n*sizeof(doublereal));

for(i=0; i<n; i++){
  for(j=0; j<n; j++){
    if(i==j){
      A[i+j*n]=4;
      if(i!=0 && i%p!=0) A[i+j*n-1]=(-1);
      if(i!=n-1 && i%p!=(p-1)) A[i+j*n+1]=(-1);
    }
    if(i==j+p) A[i+j*n]=(-1);
    if(j==i+p) A[i+j*n]=(-1);
  }
}

doublereal* CC=malloc(n*n*sizeof(doublereal));

integer pola=floor((n+1)/2);

for(i=0; i<n; i++){
    for(j=0; j<n; j++){
        if(i%2==0){
            if(j%2==0){
                CC[i/2+(j/2)*n]=A[i+j*n];
            }
            if(j%2){
                CC[i/2+(pola+(j-1)/2)*n]=A[i+j*n];
            }
        }
        if(i%2){
            if(j%2){
                CC[pola+(i-1)/2+(pola+(j-1)/2)*n]=A[i+j*n];
            }
            if(j%2==0){
                CC[pola+(i-1)/2+(j/2)*n]=A[i+j*n];
            }
        }
    }
}

doublereal* b=malloc(n*sizeof(doublereal));
doublereal* b2=malloc(n*sizeof(doublereal));
doublereal* x_T=malloc(n*sizeof(doublereal));
doublereal* x_0=malloc(n*sizeof(doublereal));
         
for(i=0; i<n; i++){
    x_T[i]=1.0;
    x_0[i]=0.0;
    b[i]=0.0;
}

char trans='N';
doublereal alpha=1.0, beta=0.0;
integer incx=1, incy=1;
dgemv_(&trans, &n, &n, &alpha, A, &n, x_T, &incx, &beta, b, &incy);
dgemv_(&trans, &n, &n, &alpha, CC, &n, x_T, &incx, &beta, b2, &incy);

doublereal aa=1;
doublereal opt=1.0;
doublereal omega;
doublereal greska, min=2.0;

for(i=0; i<10; i++){
    omega=aa+0.1*i+0.05;
    greska=sor_rjesavac (A, b, x_0, n, omega, br_iter);
    if(greska<min){
        min=greska;
        opt=omega;
    }
    printf("\n Za omega = %lf , greska je: %e .\n", omega, greska);
}
printf("\n Najbolji omega od odabranih je: %lf .\n", opt);
doublereal pom=opt;

for(i=0; i<9; i++){
    omega=pom-0.04+0.01*i;
    greska=sor_rjesavac (A, b, x_0, n, omega, br_iter);
    if(greska<min){
        min=greska;
        opt=omega;
    }
    printf("\n Za omega = %lf , greska je: %e .\n", omega, greska);
}
printf("\n Najbolji omega od odabranih je: %lf .\n", opt);
pom=opt;

for(i=0; i<9; i++){
    omega=pom-0.004+0.001*i;
    greska=sor_rjesavac (A, b, x_0, n, omega, br_iter);
    if(greska<min){
        min=greska;
        opt=omega;
    }
    printf("\n Za omega = %lf , greska je: %e \n", omega, greska);
}
printf("\n Optimalni omega je: %lf .\n", opt);

doublereal aa2=1;
doublereal opt2=1.0;
doublereal omega2;
doublereal greska2, min2=2.0;

for(i=0; i<10; i++){
    omega2=aa2+0.1*i+0.05;
    greska2=sor_rjesavac (CC, b2, x_0, n, omega2, br_iter2);
    if(greska2<min2){
        min2=greska2;
        opt2=omega2;
    }
    printf("\n Za omega = %lf , greska je: %e \n", omega2, greska2);
}
printf("\n Najbolji omega od odabranih je: %lf \n", opt2);
doublereal pom2=opt2;

for(i=0; i<9; i++){
    omega2=pom2-0.04+0.01*i;
    greska2=sor_rjesavac (CC, b2, x_0, n, omega2, br_iter2);
    if(greska2<min2){
        min2=greska2;
        opt2=omega2;
    }
    printf("\n Za omega = %lf , greska je: %e \n", omega2, greska2);
}
printf("\n Najbolji omega od odabranih je: %lf \n", opt2);
pom2=opt2;

for(i=0; i<9; i++){
    omega2=pom2-0.004+0.001*i;
    greska2=sor_rjesavac (CC, b2, x_0, n, omega2, br_iter2);
    if(greska2<min2){
        min2=greska2;
        opt2=omega2;
    }
    printf("\n Za omega = %lf , greska je: %e \n", omega2, greska2);
}
printf("\n Optimalni omega je: %lf \n", opt2);

}