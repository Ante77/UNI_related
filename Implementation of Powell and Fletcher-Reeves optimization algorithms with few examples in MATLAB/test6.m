## Funkcija Fletcher_Reeves_algorithm testira Zadatak 3 sa str. 388.
## Usput je provjereno da je (s_i)^T*H*s_j=0 za smjerove koje smo dobili.

clear
clc

f=@C8S5E3;
g=@grad_C8S5E3;
H=@Hess_C8S5E3;

epsilon=1e-5;

x01=[4;2];
x02=[3;1];
x03=[1;1];

format short e
epsilon
format short
[xmin1,it1,S1]=Fletcher_Reeves_algorithm(f,g,x01,epsilon)
[xmin2,it2,S2]=Fletcher_Reeves_algorithm(f,g,x02,epsilon)
[xmin3,it3,S3]=Fletcher_Reeves_algorithm(f,g,x03,epsilon)

## Metoda daje to?no rješenje za sve 3 po?etne iteracije za razliku od
## Powellovog osnovnog algoritma

Search_directions_are_conjugate1=S1(:,1)'*H(xmin3)*S1(:,2)
Search_directions_are_conjugate2=S2(:,1)'*H(xmin3)*S2(:,2)
Search_directions_are_conjugate3=S3(:,1)'*H(xmin3)*S3(:,2)

minimum=f(xmin1)

n=50;
F=zeros(n+1,n+1);
x1=-5;
x2=5;
y1=-5;
y2=5;
dx=(x2-x1)/n;
dy=(y2-y1)/n;
x=x1:dx:x2;
y=y1:dy:y2;

for i=1:n+1
  for j=1:n+1
    F(i,j)=f([x(i),y(j)]);
  end
end

[X,Y]=meshgrid(x,y);
figure()
surf(X,Y,F);