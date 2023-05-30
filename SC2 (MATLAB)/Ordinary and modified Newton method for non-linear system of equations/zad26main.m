clear
clc

## max 500 iteracija

F=@fja_zad2705;
df=@df_zad2705;

x0prva=[1.5;2];
x0druga=[0.5;0.4];


[x, k] = Newton_sustav (F, df, x0prva);
[x2, k2] = Newton_sustav (F, df, x0druga);
[xg, kg] = mod_newton (F, df, x0prva);
[x2g, k2g] = mod_newton (F, df, x0druga);

x
k
x2
k2
xg
kg
x2g
k2g