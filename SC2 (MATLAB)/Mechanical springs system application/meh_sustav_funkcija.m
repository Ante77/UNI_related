function f=meh_sustav_funkcija(t,y)
  t_0=0;
  a_0=9.8;
  zeta=0.1;
  omega_n=35;
  F=@(t) 0*(t<t_0)+a_0*(t>=t_0);
  
  f=[y(2);
  F(t)-2*zeta*omega_n*y(2)-omega_n^2*y(1)];
end