function y = gen_horner_trig(x,A,B)
  
  M=max(size(B));
  y=A(1)/2;
  for h=1:M-1
    y=y+A(h+1)*cos(h*x)+B(h+1)*sin(h*x);
  end
  y=y+(A(M+1)/2)*cos(M*x);
  
end