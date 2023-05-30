load ('model_portfelj_Cm.mat');

n=max(size(C));
tol=1e-8;
x0=zeros(n,1);
e=ones(n,1);

omega=sor_konvergencija(C);

[x,k,rezici]=sor(C,e,x0,tol,omega);

normae=norm(e);
[X,FLAG,RELRES,ITER,RESVEC]=pcg(C,e,1e-8,n,[],[],x0);
RESVEC=RESVEC/normae;

figure()
prvi=0:k;
drugi=0:ITER;
semilogy(drugi,RESVEC,'r',prvi,rezici,'b-')


[x2,k2,rezici2]=sor(C,m,x0,tol,omega);
normam=norm(m);
[X2,FLAG2,RELRES2,ITER2,RESVEC2]=pcg(C,m,1e-8,n,[],[],x0);
RESVEC2=RESVEC2/normam;

figure()
prvi2=0:k2;
drugi2=0:ITER2;
semilogy(drugi2,RESVEC2,'r',prvi2,rezici2,'b-')

eh=e';
meh=m';

Y=C\e;
Y2=C\m;

omegamin=(Y)/((e')*Y);
mip=0.05;

emip=(((eh*Y)*Y2-(eh*Y2)*Y)/((eh*Y)*(meh*Y2)-(eh*Y2)*(eh*Y2)))*mip+((meh*Y2)*X-(eh*Y2)*Y2)/((eh*Y)*(meh*Y2)-(eh*Y2)*(eh*Y2));

