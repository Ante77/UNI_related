## Funkcija Fletcher_Reeves_algorithm testira Primjer 8.5.6. sa str. 381.
## Ovo je ujedno rjesenje Zadatka 5 na str. 388 s provjerom da je (s_i)^T*H*s_j=0
## za smjerove koje smo dobili. 

clear
clc

f=@defquadfja;
g=@grad_defquadfja;
H=@Hess_defquadfja;

epsilon=1e-5;

x01=[2.5;2];
x02=[2;2];
x03=[1;2];

format short e
epsilon
format short
[xmin1,it1,S1]=Fletcher_Reeves_algorithm(f,g,x01,epsilon)
[xmin2,it2,S2]=Fletcher_Reeves_algorithm(f,g,x02,epsilon)
[xmin3,it3,S3]=Fletcher_Reeves_algorithm(f,g,x03,epsilon)

## Metoda daje tocno rješenje za sve 3 pocetne iteracije za razliku od
## Powellovog osnovnog algoritma

Search_directions_are_conjugate1=S1(:,2)'*H(xmin3)*S1(:,1)
Search_directions_are_conjugate2=S2(:,2)'*H(xmin3)*S2(:,1)
Search_directions_are_conjugate3=S3(:,2)'*H(xmin3)*S3(:,1)

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