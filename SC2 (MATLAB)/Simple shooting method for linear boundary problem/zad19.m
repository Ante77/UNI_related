T=@fjaT;
g=@fjag;

A=[1 0 -1; 2 3 1; 4 5 -2];
B=[1 0 0; 0 1 0; 0 0 1];
c=[3; -8; 5];

a=0;
b=2;
n=100;

[x,y,s]=odj_gadjanje_linrp(T,g,a,b,A,B,c,n);

s

figure()
plot(x,y(1,:),'-b','Linewidth',2);
hold on;
plot(x,y(2,:),'-c','Linewidth',2);
hold on;
plot(x,y(3,:),'-g','Linewidth',2);
legend('y_1(x)','y_2(x)','y_3(x)');
xlabel('x');
ylabel('y_i(x)');
hold off;