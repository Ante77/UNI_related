function [x,y] = odj_pred_kor (f, a, b, y0, n, m)

h=(b-a)/n;
bc=a+3*h;
[x,y]=odj_rk4(f,a,bc,y0,3);

F=zeros(n+1,2);
F(1:4,1)=x;
F(1:4,2)=y;
  
for i=4:n
  F(i+1,2)=F(i,2)+h/24*(55*f(F(i,1),F(i,2))-59*f(F(i-1,1),F(i-1,2))+37*f(F(i-2,1),F(i-2,2))-9*f(F(i-3,1),F(i-3,2)));
  F(i+1,1)=F(i,1)+h;
  
  for j=1:m
    F(i+1,2)=F(i,2)+h/24*(9*f(F(i+1,1),F(i+1,2))+19*f(F(i,1),F(i,2))-5*f(F(i-1,1),F(i-1,2))+f(F(i-2,1),F(i-2,2)));
  endfor
endfor
  
x=F(:,1);
y=F(:,2);
  
endfunction