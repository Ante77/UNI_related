load('K.txt');

M=zeros(4,4);
M(1,1)=-1/sqrt(2);
M(2,2)=-1/sqrt(5);
M(3,3)=-1/sqrt(3);
M(4,4)=-1/sqrt(6);
x0=[0.5;0.5;0.5;0.5];

%Ap=(M^(1/2))\K/(M^(1/2));
Ap=M*K*M;
% nema minres funkcije u Octavi
I=eye(4);
k=0;
b=x0;

A=Ap-4*I;
x=1:4;

maxit=20;
tol=1e-8;

figure()
u=[0.5;0.5;0.5;0.5];
ro=((u')*Ap*u)/((u')*u);
rezidual(1)=norm(Ap*u-ro*u);
i=1;

while rezidual(i)>1e-8
    i=i+1;
    b=u;
	[u,flag,relres,iter,resvec]=minres(A,b,tol,maxit,[],[],x0);
	k=k+1;
	axis([0.5 4.5 -1 1])
	xlabel('Indeks komponente');
	ylabel('Komponenta');
	title(sprintf('Aproksimacija svojstvenog vektora: iteracija %d',k));
	set(gca,'XTick',1:4);
	set(gca,'XTickLabel',{'1','2','3','4'});
    
    u=u/norm(u);
    ro=((u')*Ap*u);
    rezidual(i)=norm(Ap*u-ro*u);
    
    plot(x,u,'b-');
	grid on;
	pause(0.5);
end

figure()
p=0:4;
semilogy(p,rezidual,'m*-');
xlabel('Broj iteracije');
ylabel('Norma reziduala')

[V,D]=eig(Ap);