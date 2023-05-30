h=0.05;
x=0:0.05:1;
y=0:0.05:1;

A=gallery('poisson',19);
A=(1/(h*h))*A;
b=zeros(361,1);
b(181)=10000;

tol=1e-8;
maxit=500;
u0=zeros(361,1);
bnorm=norm(b);

c1=cond(full(A));
[u,flag,relres,iter,resvec]=pcg(A,b,tol,maxit,[],[],u0);
resvec=resvec/bnorm;

R=cholinc(A,'0');
c2=cond(full(R'\A/R));
[u2,flag2,relres2,iter2,resvec2]=pcg(A,b,tol,maxit,R,R',u0);
resvec2=resvec2/bnorm;

figure()
prvi=0:iter;
drugi=0:iter2;
semilogy(prvi,resvec,'r',drugi,resvec2,'b-')
legend('obicna','prekondicionirana');

U=reshape(u,19,19);
figure()
x=0.05:0.05:0.95;
y=0.05:0.05:0.95;
surf(x,y,U);