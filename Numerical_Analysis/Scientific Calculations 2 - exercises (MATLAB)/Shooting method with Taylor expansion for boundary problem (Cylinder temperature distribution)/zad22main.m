f=@odj_primjer_distr_temp_f3;
n=100;
a=0;
b=1;
alfa=0;
beta=0;
s0=1;

[x,y,s,k]=odj_gadjanje_2rru(f,a,b,alfa,beta,s0,n);

s
k

figure()
plot(x,y(1,:),'-b');
legend('aproksimativno rjesenje');
xlabel('x');
ylabel('y');