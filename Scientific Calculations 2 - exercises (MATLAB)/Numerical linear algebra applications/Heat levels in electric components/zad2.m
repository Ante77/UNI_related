n=7;
tol=1e-8;
x0=ones(n,1);

A=[1,0,0,0,-0.5,0.5,0;
    0,1,0,0,0,-2,2;
    0,0,1,0,-1/35,0,0;
    0,0,0,1,0,-1/0.7,0;
    0,1,0,0,0,0,-1;
    1,0,1,0,0,0,0;
    1,-1,0,-1,0,0,0];
b=[0;0;-20/35;-20/0.7;-20;10;0];

[X,FLAG,RELRES,ITER,RESVEC]=gmres(@mdAx,b,[],tol,n,[],[],x0);

normab=norm(b);
RESVEC=RESVEC/normab;

RESVEC(8)=norm(A*X-b)/normab;

p=0:7;

semilogy(p,RESVEC,'g-');
xlabel('Broj iteracije');
ylabel('Relativna norma reziduala');
grid;

% u iducem zadatku matricu sustava podijelit s h^2  