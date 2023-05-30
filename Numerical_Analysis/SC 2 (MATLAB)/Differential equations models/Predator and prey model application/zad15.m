f=@fja_zad_15;

p0=[5000;100];
n=300;
tol=1e-2;
a=0;
b=30;

[x,y]=odj_rk23v(f,a,b,p0,tol);
[x1,y1]=odj_rk4v(f,a,b,p0,n);

figure()

subplot(2,1,1);
plot(x,y(1,:),'-ro');
hold on;
plot(x1,y1(1,:),'-b+');
legend('rk23','rk4');
title('Populacija plijena');
xlabel('Vrijeme');
ylabel('Broj plijena');
hold off;

subplot(2,1,2);
plot(x,y(2,:),'-ro');
hold on;
plot(x1,y1(2,:),'-b+');
legend('rk23','rk4');
title('Populacija grabezljivca');
xlabel('Vrijeme');
ylabel('Broj grabezljivca');
hold off;