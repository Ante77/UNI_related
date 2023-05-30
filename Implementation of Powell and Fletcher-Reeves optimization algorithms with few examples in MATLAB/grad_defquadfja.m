function [rez]= grad_defquadfja (v)
  rez=[2*(2*v(1)-v(2))*2;
  -2*(2*v(1)-v(2))+2*(v(2)+1)];
end