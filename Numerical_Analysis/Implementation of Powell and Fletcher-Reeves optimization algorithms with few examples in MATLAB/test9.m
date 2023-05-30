## Usporedba kako Quadratic Search radi brze od Matlab funkcije "fminbnd" ako jet
## potrebna velika tocnost.ž

## Powellov osnovni algoritam na Rosenbrockovoj funkciji 4 varijable
## LOSE RADI!

clear
clc

f=@Rosenbrock_fja;

epsilon1=1e-2;
epsilon2=1e-5;
epsilon3=1e-10
x01=[0;0;0;0];
x02=[1;2;3;4];

## Ne konvergiraju sve k istim vrijednostima... (1,1,1,1) je globalni minimum, ne znam je li jedini lokalni.
## Neke konvergiraju k slicnima

format short
## Quadratic Search
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x01,epsilon1) # 10 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x01,epsilon2) # 25 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x01,epsilon3) # 27 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x02,epsilon1) # 12 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x02,epsilon2) # 17 iter
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x02,epsilon3) # 46 iter
## fminbnd
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x01,epsilon1) # 11 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x01,epsilon2) # 13 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x01,epsilon3) # 334 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x02,epsilon1) # 2 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x02,epsilon2) # 14 iter
[xmin2,it2]=Powell_algorithm2_bez_ispisa(f,x02,epsilon3) # 1802 iter