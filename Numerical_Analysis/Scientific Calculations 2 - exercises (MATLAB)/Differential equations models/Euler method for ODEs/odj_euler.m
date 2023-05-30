function [x,y] = odj_euler (f, a, b, y0, n)

h=(b-a)/n;
x=zeros(n+1,1);
y=zeros(n+1,1);
x(1)=a;
y(1)=y0;

for i=1:n
    x(i+1)=a+i*h;
    y(i+1)=y(i)+h*f(x(i),y(i));
end

end