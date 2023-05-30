function f=term_obrada_funcija(t,T)
  V=(0.01/2)^2*pi;
  ro=7822;
  m=ro*V;
  c=444;
  eps=0.7;
  H_1=15;
  H_2=100;
  T_0=294.15;
  t_h=70;
  Q_e=@(t) 3000*(t<t_h);
  H=@(t) H_1*(t<t_h)+H_2*(t>=t_h);
  A_p=2*(0.01/2)*pi;
  sigma=5.67*(1e-8);
  
  f=(1/(m*c))*(Q_e(t)-A_p*(H(t)*(T-T_0)+eps*sigma*(T^4-T_0^4)));
end