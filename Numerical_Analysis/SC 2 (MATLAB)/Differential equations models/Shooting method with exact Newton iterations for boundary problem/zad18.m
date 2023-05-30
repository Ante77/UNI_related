f=@fja1;
dfy=@df2;
dfdy=@dfdy2;
n=100;
a=0;
b=1;
alfa=4;
beta=1;
s01=-4;
s02=-20;

[x,y,s,k]=odj_gadjanje_2rru_en(f,dfy,dfdy,a,b,alfa,beta,s01,n);
[x1,y1,s1,k1]=odj_gadjanje_2rru_en(f,dfy,dfdy,a,b,alfa,beta,s02,n);

s
k
s1
k1

figure()
plot(x,y(1,:),'-b');
hold on;
plot(x1,y1(1,:),'-r');
hold off;