f=@term_obrada_funkcija;

a=0;
b=210;
tol=1e-3;
y0=294.15;
n=33;

[x,y]=odj_rk23 (f,a,b,y0,tol);
[x1,y1]=odj_rk4 (f,a,b,y0,n);

figure()
plot(x,y,'-or');
hold on;
plot(x1,y1,'-+b');
xlabel('Vrijeme (s)');
ylabel('Temperatura (K)');
legend('rk23','rk4');
hold off;

[x2,y2]=odj_rk4 (f,a,b,y0,2*n);

figure()
plot(x,y,'-or');
hold on;
plot(x2,y2,'-+b');
xlabel('Vrijeme (s)');
ylabel('Temperatura (K)');
legend('rk23','rk4 h/2');
hold off;