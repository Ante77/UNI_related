#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"
#include "GMRES.c"

integer size=100;
doublereal* A;

int matvec (doublereal* alpha, doublereal* x, doublereal* beta, doublereal* y){
    char trans='N';
    integer incx=1, incy=1;

    dgemv_(&trans, &size, &size, alpha, A, &size, x, &incx, beta, y, &incy);
    return 0;
}

int psolve (doublereal* x, doublereal* b){
    integer incx=1, incb=1;
    
    dcopy_(&size, b, &incx, x, &incb);
    return 0;
}

int main(integer argc, char*argv[]){

    int brojac=1;

    A=malloc(size*size*sizeof(doublereal));

    integer iseed[]={401,41,777,69};
    integer idist=1;
    integer nxn=size*size;

    doublereal *V=malloc(size*size*sizeof(doublereal));

    dlarnv_(&idist, iseed, &nxn, V);

    doublereal *tau=malloc(size*size*sizeof(doublereal));
    doublereal *work=malloc(size*size*sizeof(doublereal));
    integer info;
    integer broj_hhr=size-1;

    dgeqrf_(&size, &size, V, &size, tau, work, &size, &info);
    dorgqr_(&size, &size, &size, V, &size, tau, work, &nxn, &info);

    integer size_of_f=size+1;

    doublereal *f=malloc(size_of_f*sizeof(doublereal));
    doublereal *g=malloc(size*sizeof(doublereal));
    doublereal *b=malloc(size*sizeof(doublereal));

    integer i, j;

    for(i=0; i<101; i++)
        f[i]=100-i;
    /*
    for(i=0; i<100; i++)
        g[i]=sqrt(f[i]*f[i]-f[i+1]*f[i+1]);
*/

    for(i=0; i<100; i++)
        g[i]=sqrt((100-i-1+1)*(100-i-1+1)-(100-i-1)*(100-i-1));

/*
    printf("Dokle smo dosli?? %d \n", brojac++);
    for(i=0; i<=100; i++){
        printf(" %lf ", f[i]);
    }
    printf("\n");
    */
 /*   printf("Dokle smo dosli?? %d \n", brojac++);
    for(i=0; i<100; i++){
        printf(" %lf ", g[i]);
    }
*/
    printf("\n");
    doublereal alpha=1.0, beta=0.0;
    char trans='N';
    integer incx=1, incy=1;
    for(j=0; j<100; j++)
            b[j]=0.0;
    dgemv_(&trans, &size, &size, &alpha, V, &size, g, &incx, &beta, b, &incy);

    doublereal *B=malloc(size*size*sizeof(doublereal));

    for(i=0; i<100; i++){
        for(j=0; j<100; j++)
            A[i+j*size]=0.0;
    }

/*      zlu ne tribalo...
    char uplo='A';
    dlacpy_(&uplo, &size, &size, V, &size, B+100, &size);
*/
    printf("\n \n OVDJE b \n \n");
    for(i=0; i<100; i++){
        B[i]=b[i];
        printf(" %lf ", b[i]);
    }



    printf("\n Norma od b: %lf \n", dnrm2_(&size, b, &incx));
    for(i=100; i<100*100; i++){
        B[i]=V[i-100];
    }


    doublereal *pom=malloc(size*size*sizeof(doublereal));

    for(i=0; i<10000; i++){
        pom[i]=0.0;
    }

    doublereal *pom2=malloc(size*size*sizeof(doublereal));

    for(i=0; i<10000; i++){
        pom2[i]=0.0;
    }

    for(i=0; i<100; i++){
        pom[i+i*size+1]=1.0;
    }
    pom[9900]=1.0;

    printf(" \n OVDJE B matrica \n \n ");
// mnozim B*ono cudo prvo
 /*   for(i=0; i<100; i++){
        for(j=0; j<100; j++){
            printf(" %lf ", B[i+j*size]);
        }
        printf("\n");
    }
 */ 
    integer transa='N', transb='N';
    doublereal alfa=1.0, betha=0.0;
    dgemm_(&transa, &transb, &size, &size, &size, &alfa, B, &size, pom, &size, &betha, pom2, &size);
// pokusaj inverza od B

    doublereal *work2=malloc(size*size*sizeof(doublereal));
    integer *ipiv=malloc(size*sizeof(integer));
    integer info2, info3;
    dgetrf_(&size, &size, B, &size, ipiv, &info2);

    dgetri_(&size, B, &size, ipiv, work2, &size, &info3);

// mnozim umnozak prve dvi matrice s inverzom od B    
    dgemm_(&transa, &transb, &size, &size, &size, &alfa, pom2, &size, B, &size, &betha, A, &size);

    doublereal *x_0=malloc(size*sizeof(doublereal));
    for(i=0; i<100; i++)
        x_0[i]=0.0;
    integer restrt=size;
    doublereal *work3=malloc(size*(restrt+4)*sizeof(doublereal));
    integer ldw=100;
    doublereal *h=malloc((size+1)*(restrt+2)*sizeof(doublereal));
    integer ldh=101;
    integer iter=100;
    doublereal resid=1e-5;
    integer infor;

    printf("\n \n OVDJE b \n \n");
    for(i=0; i<100; i++){
        printf(" %lf ", b[i]);
    }
    printf("\n \n OVDJE x_0 \n \n");
    for(i=0; i<100; i++){
        printf(" %lf ", x_0[i]);
    }

    printf(" \n OVDJE A matrica \n \n ");
// mnozim B*ono cudo prvo
  /*  for(i=0; i<100; i++){
        for(j=0; j<100; j++){
            printf(" %lf ", pom[i+j*size]);
        }
        printf("\n");
    }
    printf("\n \n");
*/

    gmres_(&size, b, x_0, &restrt, work3, &ldw, h, &ldh, &iter, &resid, matvec, psolve, &infor);

    printf("Afrika paprika?");

    integer lwork1=3*size;
    double* wr=malloc(size*sizeof(double));//out realni dio
    double* wi=malloc(size*sizeof(double));//out imaginarni dio
    double* vl=malloc(size*sizeof(double)); 
    double* vr=malloc(size*sizeof(double));
    double* work1=malloc(lwork1*sizeof(double));
    
    dgeev_(&transa,&transa,&size,A,&size,wr,wi,vl,&size,vr,&size,work1,&lwork1,&info);

    double* aaa=malloc(size*sizeof(double));
    double* bbb=malloc(size*sizeof(double));
    double* stotar=malloc(size*sizeof(double));
    double* stotai=malloc(size*sizeof(double));
    double* stota=malloc(size*sizeof(double));
/*    for(i=0; i<100; i++){
        aaa[i]=sqrt(wr[i]*wr[i]+wi[i]*wi[i]);
        aaa[i]=pow(aaa[i],100);
        printf("\n %lf \n",aaa[i]);
    }
*/
    for(i=0; i<size; i++){
        aaa[i]=atan2(wr[i],wi[i]);
        doublereal stota;
        stotar[i]=cos(100*aaa[i]);
        stotai[i]=sin(100*aaa[i]);
        stota=sqrt(stotar[i]*stotar[i]+stotai[i]*stotai[i]);

        printf("\n %lf  Re: %lf  Im: %lf \n", stota, stotar[i], stotai[i]);
    }

        
    return 0;

}

