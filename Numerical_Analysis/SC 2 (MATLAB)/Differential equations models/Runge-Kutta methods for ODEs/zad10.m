f = @(x, y) -y-5*exp(-x)*sin(5*x);

y0=1;
[x,y]=odj_rk4(f,0,3,y0,30);

figure()
plot(x,y,'*--');
xlabel('x');
ylabel('y');
hold on;
fplot('exp(-x)*cos(5*x)',[0,3],'r-');
hold off;
legend('aproksimativno','egzaktno');

rez=zeros(31,1);
for i=1:31
    rez(i)=abs(y(i)-exp(-0.1*(i-1))*cos(5*(0.1*(i-1))));
end
m=max(rez);