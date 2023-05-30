function [opt] = sor_konvergencija(A)

D=diag(diag(A));
R=triu(A,1);
L=tril(A,-1);
min=1;
opt=1;
br=1;
for omega=0:0.01:2
    T=(inv(D+omega*L))*((1-omega)*D-omega*R);
    sprovi(br)=max(abs(eig(T)));
    if sprovi(br)<min
        opt=omega;
        min=sprovi(br);
    end
    br=br+1;
end
omega=0:0.01:2;
plot(omega,sprovi,'r-');
axis( [0 2 0 1]);
xlabel('omega');
ylabel('spr (T_omega)');

end