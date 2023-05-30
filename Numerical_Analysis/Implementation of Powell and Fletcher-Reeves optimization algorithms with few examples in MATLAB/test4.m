## Funkcije Powell_algorithm_bez_ispisa i Powell_algorithm2 testiraju 
## Zadatak 3 sa str. 388.

clear
clc

f=@C8S5E3;

epsilon=1e-5;

x01=[4;2];
x02=[3;1];
x03=[1;1];

format short e
epsilon
format short
[xmin1,it1]=Powell_algorithm_bez_ispisa (f,x01,epsilon)
[xmin2,it2]=Powell_algorithm_bez_ispisa (f,x02,epsilon)
[xmin3,it3]=Powell_algorithm_bez_ispisa(f,x03,epsilon)

[xmin1mat,it1mat]=Powell_algorithm2 (f,x01,epsilon)
[xmin2mat,it2mat]=Powell_algorithm2 (f,x02,epsilon)
[xmin3mat,it3mat]=Powell_algorithm2 (f,x03,epsilon)

## Iz po?etne to?ke (1,1) ?emo nakon prvog prolaska kroz while petlju dobiti 2 kolinearna smjera
## s_1 i s_2 što ?e rezultirati s vrijednosti alfa_2=0 jer smo ve? prije pronalaska smjera s_2
## bili na minimumu funkcije u pravcu toga smjera. Zna?i da se sljede?a iteracija ne?e promijeniti
## i while petlja ?e se završiti. 

## Taj "problem" se može riješiti Fletcher-Reeves metodom.

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