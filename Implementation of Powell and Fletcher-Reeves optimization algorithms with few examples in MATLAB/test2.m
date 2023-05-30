## Funkcija Powell_algorithm_bez_ispisa testira Primjer 8.5.6. sa str. 381.
## Za minimiziranje funkcija f(alfa)=f(x_i+alfa*s_i) i f(alfa)=f(x_n+alfa*s_n) koristi
## funkciju Quadratic_search_algorithm koja implementira Quadratic search metodu za
## minimizaciju funkcije jedne varijable.

## Funkcija ispisuje samo to?ke minimuma x i broj prolazaka kroz while petlju, tj. broj
## iteracija. Tolerancija epsilon ne utje?e zna?ajno na rješenje jer se u danom primjeru
## skoro uvijek pogodi to?an smjer.

## Pošto se uvjet za "konvergenciju" da je razlika izme?u dvije uzastopne iteracije
## manja od epsilon, uvijek se izvrši jedan korak petlje više nego što je napisano u knjizi.

clear
clc

f=@defquadfja;

epsilon1=1e-1;
epsilon2=1e-2;
epsilon3=1e-5;
epsilon4=1e-10;
x01=[2.5;2];
x02=[2;2];
x03=[1;2];

format short e
epsilon1
format short
[xmin1,it1]=Powell_algorithm_bez_ispisa (f,x01,epsilon1)
[xmin2,it2]=Powell_algorithm_bez_ispisa (f,x02,epsilon1)
[xmin3,it3]=Powell_algorithm_bez_ispisa (f,x03,epsilon1)

format short e
epsilon2
format short
[xmin1,it1]=Powell_algorithm_bez_ispisa (f,x01,epsilon2)
[xmin2,it2]=Powell_algorithm_bez_ispisa (f,x02,epsilon2)
[xmin3,it3]=Powell_algorithm_bez_ispisa (f,x03,epsilon2)

format short e
epsilon3
format short
[xmin1,it1]=Powell_algorithm_bez_ispisa (f,x01,epsilon3)
[xmin2,it2]=Powell_algorithm_bez_ispisa (f,x02,epsilon3)
[xmin3,it3]=Powell_algorithm_bez_ispisa (f,x03,epsilon3)

format short e
epsilon4
format short
[xmin1,it1]=Powell_algorithm_bez_ispisa (f,x01,epsilon4)
[xmin2,it2]=Powell_algorithm_bez_ispisa (f,x02,epsilon4)
[xmin3,it3]=Powell_algorithm_bez_ispisa (f,x03,epsilon4)

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