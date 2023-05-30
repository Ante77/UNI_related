load('model_orbita_polozaji.mat');

n=5;
A=zeros(5,5);
A=praviredak(A,p1,1);
A=praviredak(A,p2,2);
A=praviredak(A,p3,3);
A=praviredak(A,p4,4);
A=praviredak(A,p5,5);
b=zeros(n,1);
b(1)=-p1(1)*p1(1);
b(2)=-p2(1)*p2(1);
b(3)=-p3(1)*p3(1);
b(4)=-p4(1)*p4(1);
b(5)=-p5(1)*p5(1);

a=A\b;

figure()
x=-5:0.1:5;
y=x;
[X,Y]=meshgrid(x,y);
Z=fja(X,Y,a);
contour(x,y,Z,[0 0],'r');
grid on;

figure()
surf(x,y,Z);
grid on;