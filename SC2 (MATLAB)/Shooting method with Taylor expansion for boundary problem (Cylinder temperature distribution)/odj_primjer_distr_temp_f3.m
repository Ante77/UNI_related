function f=odj_primjer_distr_temp_f3(x,y)
  alfa=0.8;
  if x>=0 && x<=1e-2
    f=[y(2);
    1/2*alfa*exp(y(3))*(-1+3/8*x^2*alfa*exp(y(3)));
    0];
  endif
  
  if x>=1e-2 && x<=1
    f=[y(2);
    -y(2)/x-alfa*exp(y(1));
    0];
  endif
end