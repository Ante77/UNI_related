load('model_dekonvolucija_uy.mat');

n=18;
m=102;
U=zeros(m+1,n+1);

% u(i)=up(i+n+1)
% prvi red u ide od 0 do -n, tj. od n+1 do 1 po up

for i=1:m+1
    for j=1:n+1
        U(i,j)=up(i+1-j+n);
    end
end

% y(i)=yp(i+1)
y=zeros(m+1,1);

for i=1:m+1
    y(i)=yp(i);
end

[u,s,v]=svd([U y]);
[u2,s2,v2]=svd([U]);
sigma_n=s2(19,19)
sigma_n_plus_jedan=s(20,20)

load('model_dekonvolucija_hh.mat');

T1=-1*eye(19);
V12=v(1:19,20);
V22=v(20,20);
T2=eye(1);

h=T1*V12/V22/T2;

figure()
p=1:19;
plot(p,hh,'g-',p,h,'r-');
legend('egzaktno','aproksimativno')

% ZA SLJEDECI PUT:
% f=@(x) sin(5*x*x)
% f=@(x,y) -y-5*exp(-x).*sin(5*x*x)