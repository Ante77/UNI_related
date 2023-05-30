function [min,it,Sovi]=Fletcher_Reeves_algorithm(f,g,x0,epsilon)
  dim=max(size(x0));
  s=zeros(dim,1);
  beta=1;
  G0=g(x0);
  it=0;
  x=x0;
  G1=G0;
  
  while norm(G1)>=epsilon
    s=-G1+beta*s;
    Sovi(:,it+1)=s;
    fja=@(alfa) f(x+alfa*s);
    alfamin=fminbnd(fja,-100,100);
    x=x+alfamin*s;
    G1=g(x);
    beta=(norm(G1)^2)/(norm(G0)^2);
    G0=G1;
    it=it+1;
    norm(G1);
  end
  min=x;
end
