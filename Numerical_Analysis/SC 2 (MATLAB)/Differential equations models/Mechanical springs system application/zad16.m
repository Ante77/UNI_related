f=@meh_sustav_funkcija;

y0=[0;0];
a=0;
b=1.5;
tol=1e-5;

[x,y]=odj_rk23v(f,a,b,y0,tol);

figure()

subplot(2,1,1);
plot(x,y(1,:),'-ro');
legend('rk23');
xlabel('Vrijeme');
ylabel('Polozaj');
title('Graf polozaja');

subplot(2,1,2);
plot(x,y(2,:),'-ro');
legend('rk23');
xlabel('Vrijeme');
ylabel('Brzina');
title('Graf brzine');