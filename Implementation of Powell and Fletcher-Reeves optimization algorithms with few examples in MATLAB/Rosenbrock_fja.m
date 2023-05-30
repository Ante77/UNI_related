function [rez]=Rosenbrock_fja(v)
  rez=100*(v(1)^2-v(2))^2+(1-v(1))^2+100*(v(3)^2-v(4))^2+(1-v(3))^2;
end