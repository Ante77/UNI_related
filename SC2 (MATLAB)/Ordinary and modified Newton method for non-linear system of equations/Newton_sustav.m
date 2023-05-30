function  [x, k] = Newton_sustav (f, df, x0)
  
  dim=max(size(x0));
  N=500;
  k=1;
  F=zeros(dim,1);
  J=zeros(dim,dim);
  F=f(x0);
  J=df(x0);
  
  kor=J\F;
  x=x0-kor;
  
  flag=0;
  
  if (norm(x-x0)/norm(x0))<eps
    flag=1;
  end
  
  if flag==0
    while (norm(x-x0)/norm(x0))>=eps && k<N
      x0=x;
      F=f(x);
      J=df(x);
      kor=J\F;
      x=x-kor;
      k=k+1;
    end
  end

end