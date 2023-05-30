load('model_gravitacija_g.mat');

%fplot('sin(pi*x)+0.5*sin(2*pi*x)',[0,1],'r-');

a=0;
b=1;
d=0.25;
m=15;
n=15;
h=(b-a)/(n-1);
x=0:h:1;


k=zeros(m,n);
K=zeros(m,n);
s=0:h:1;
t=0:h:1;

for i=1:m
    for j=1:n
        k(i,j)=d/((d^2+(s(i)-t(j))^2)^(3/2));
    end
end

for i=1:m
    for j=2:n-1
        K(i,j)=((b-a)/(n-1))*k(i,j);
    end
    K(i,1)=((b-a)/(2*(n-1)))*k(i,1);
    K(i,n)=((b-a)/(2*(n-1)))*k(i,n);
end

[U,S,V]=svd(K);

s_vrij=diag(S);

figure()
plot(s_vrij,'bo')

fun=zeros(m,1);
for i=1:15
    fun(i)=sin(pi*x(i))+0.5*sin(2*pi*x(i));
end
figure()
plot(x,fun,'b-');
axis([0 1 0 1.4])

greska=zeros(m,1);

for p=1:15
    xtilda=V(:,1:p)*(S(1:p,1:p)\((U(:,1:p)')*g));
    plot(x,xtilda,'r--');
    hold on;
    plot(x,fun,'b-');
    xlabel('x');
    ylabel('y');
    title(sprintf('Aproksimacija ranga %d',p));
    axis([0 1 0 1.4])
    pause(0.5);
    hold off;
    greska(p)=norm(xtilda-fun);
end

[razlika,indeks]=min(greska)