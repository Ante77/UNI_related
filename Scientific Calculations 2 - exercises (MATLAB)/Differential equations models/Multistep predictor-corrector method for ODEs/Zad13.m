f = @(x,y) -y-5.*exp(-x)*sin(5.*x);
egz = @(x) exp(-x).*cos(5.*x);
a=0;
b=3;
y0=1;

n=30;
m=3;

[x30,y30]=odj_pred_kor(f,a,b,y0,n,m);

figure()
plot(x30,egz(x30),'r-');
hold on;
plot(x30,y30,'b-*');
hold off;
legend('egzaktno rjesenje','prediktor-korektor');
title('n=30 , m=3');

n=300;
m=3;

[x300,y300]=odj_pred_kor(f,a,b,y0,n,m);

figure()
plot(x300,egz(x300),'r-');
hold on;
plot(x300,y300,'b-*');
hold off;
legend('egzaktno rjesenje','prediktor-korektor');
title('n=300 , m=3');

error1_n30_m3=max(abs(egz(x30)-y30))
error2_n300_m3=max(abs(egz(x300)-y300))