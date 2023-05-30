function [x,y,s]=odj_gadjanje_linrp(T,g,a,b,A,B,c,n)
% Y'=T(x)*Y, Y1'=[T11 T12 T13]*[Y11;Y21;Y31];
f=@(x,y) T(x)*y;
%[T(x)(1,1)*x(1)+T(x)(1,2)*x(2)+T(x)(1,3)*x(3);
 %   T(x)(2,1)*x(1)+T(x)(2,2)*x(2)+T(x)(2,3)*x(3);
  %  T(x)(3,1)*x(1)+T(x)(3,2)*x(2)+T(x)(3,3)*x(3)];
y01=[1;0;0];
y02=[0;1;0];
y03=[0;0;1];

[x1,y1]=odj_rk4v(f,a,b,y01,n);
[x2,y2]=odj_rk4v(f,a,b,y02,n);
[x3,y3]=odj_rk4v(f,a,b,y03,n);

Y(:,1)=y1(:,n+1);
Y(:,2)=y2(:,n+1);
Y(:,3)=y3(:,n+1);

DF0=A+B*Y;

f2=@(x,y) T(x)*y+g(x);
[x4,y4]=odj_rk4v(f2,a,b,[0;0;0],n);

F0=B*y4(:,n+1)-c;

s=-(DF0\F0);

[x,y]=odj_rk4v(f2,a,b,s,n);

end