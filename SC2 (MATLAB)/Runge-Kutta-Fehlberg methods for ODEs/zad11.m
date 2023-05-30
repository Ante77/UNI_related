f = @(x,y) (1-y)*(x>=0 & x<=pi)-5*y*(x>pi & x<=4);
egz = @(x) (x>=0 & x<=pi).*(1-exp(-x))+(x>pi & x<=4).*((1-exp(-pi))*exp(-5*(x-pi)));

y0=0;
tol=1e-5;

[x2,y2]=odj_rk2(f,0,4,y0,37);
[x23,y23]=odj_rk23(f,0,4,y0,tol);

plot(x23,y23,'r*-');
hold on;
plot(x2,y2,'bo-');
hold on;
fplot(egz,[0,4]);
hold off;
xlabel('x');
ylabel('y');
legend('RK23','RK2','egzaktno');

cv23=egz(x23(2:38));
cv2=egz(x2(2:38)');

errrk23=max(abs(cv23-y23(2:38)))
errrk22=max(abs(cv2-y2(2:38)'))