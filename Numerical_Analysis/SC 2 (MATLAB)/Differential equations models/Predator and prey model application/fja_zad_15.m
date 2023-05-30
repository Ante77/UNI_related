function f=fja_zad_15(t,p)
  alfa_1=2;
  delta_1=0.02;
  alfa_2=0.0002;
  delta_2=0.8;
  
  f=[alfa_1*p(1)-delta_1*p(1)*p(2);
  alfa_2*p(1)*p(2)-delta_2*p(2)];
end