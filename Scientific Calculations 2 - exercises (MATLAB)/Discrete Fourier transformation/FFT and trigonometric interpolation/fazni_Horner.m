function y = fazni_Horner(x,bete)
  
  epsilon=exp(x*sqrt(-1));
  y=bete(max(size(bete)));
  for j=max(size(bete))-1:-1:1
    y=y*epsilon+bete(j);
  end
  
end