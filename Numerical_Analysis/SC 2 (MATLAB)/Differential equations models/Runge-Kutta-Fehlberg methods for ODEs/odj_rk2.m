function [x,y] = odj_rk2 (f, a, b, y0, n)

h=(b-a)/n;
x=zeros(n+1,1);
y=zeros(n+1,1);
x(1)=a;
y(1)=y0;
omega1=214/891;
omega2=1/33;
omega3=650/891;

for i=1:n
    x(i+1)=a+i*h;
    k1=f(x(i),y(i));
    k2=f(x(i)+1/4*h,y(i)+h*1/4*k1);
    k3=f(x(i)+27/40*h,y(i)+h*((-189/800)*k1+729/800*k2));
    y(i+1)=y(i)+h*(omega1*k1+omega2*k2+omega3*k3);
endfor

endfunction