## Usporedba kako Quadratic Search radi brze od Matlab funkcije "fminbnd" ako jet
## potrebna velika tocnost.

## Powellov osnovni algoritam na Woodovoj funkciji 4 varijable
## LOSE RADI!

clear
clc

f=@Wood_function;

epsilon1=1e-2;
epsilon2=1e-5;
epsilon3=1e-10
x01=[0;0;0;0];
x02=[1;2;3;4];

## Ne konvergiraju sve k istim vrijednostima... (1,1,1,1) je globalni minimum, ne znam je li jedini lokalni.
## Neke konvergiraju k slicnima

format short
## Quadratic Search
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x01,epsilon1) # 4 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x01,epsilon2) # 11 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x01,epsilon3) # 14 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x02,epsilon1) # 9 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x02,epsilon2) # 11 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x02,epsilon3) # 11 iter
## fminbnd
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x01,epsilon1) # 4 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x01,epsilon2) # 9 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x01,epsilon3) # 67 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x02,epsilon1) # 4 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x02,epsilon2) # 131 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x02,epsilon3) # 343 iter