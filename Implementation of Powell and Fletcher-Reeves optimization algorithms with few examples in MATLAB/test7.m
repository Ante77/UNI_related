## Woodova funkcija 4 varijable

## Funkcija Fletcher_Reeves_algorithm2 u svakoj iteraciji ispisuje gradijent(x) koji je
## uvjet konvergencije.
## Vrijednost grad(x) se zanimljivo ponaša u ovisnosti u po?etnoj iteraciji.
## Norma gradijenta se ?esto po?ne spuštati i onda eksplodira opet.

## EGZAKTNO RJEŠENJE JE [1;1;1;1] (globalni minimum)

clear
clc

f=@Wood_function;
g=@grad_Wood_function;

epsilon=1e-3;

## Iako metoda ne konvergira daljih pocetnih iteracija, ponasanje grad(x) koje se ispisuje
## je zanimljivo.

x01=[0;0;0;0]; ## ne konvergira u razumnom vremenu
x02=[1;2;3;4]; ## ne konvergira u razumnom vremenu
x03=[1;1.5;0.5;1]; ## 99 iter
x04=[0.5;1;1;1.5]; ## ne konvergira u razumnom vremenu
x05=[0.8;1;1;1.2]; ## 18 iter

format short e
epsilon
format short
[xminP,itP]=Powell_algorithm_bez_ispisa(f,x01,epsilon)
[xminF03,itF03]=Fletcher_Reeves_algorithm2(f,g,x03,epsilon)
[xminF05,itF05]=Fletcher_Reeves_algorithm2(f,g,x05,epsilon)
##[xminF04,itF04]=Fletcher_Reeves_algorithm2(f,g,x04,epsilon)
##[xminF02,itF02]=Fletcher_Reeves_algorithm2(f,g,x02,epsilon);
##[xminF,itF]=Fletcher_Reeves_algorithm2(f,g,x01,epsilon);