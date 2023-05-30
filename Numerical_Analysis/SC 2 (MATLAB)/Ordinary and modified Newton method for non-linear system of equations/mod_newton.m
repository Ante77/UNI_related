function  [x, k] = mod_newton (f, df, x0)

  dim=max(size(x0));
  N=500;
  k=1;
  F=zeros(dim,1);
  J=zeros(dim,dim);
  F=f(x0);
  J=df(x0);
  
  x=x0;
  d=J\F;
  gama=1/(cond(J));
  
  h=@(x) f(x)'*f(x);
  hk=@(tau) h(x-tau*d);
  
    nasao=0;
    j=0;
    while(nasao==0)
        lijevi=hk(2^(-j));
        desni=hk(0)-(2^(-j))*(gama/4)*norm(d)*norm(2*F'*J);
        if lijevi<=desni
            nasao=1;
        end
        j=j+1;
    end
    min=hk(1);
    mini=0;
    for i=0:j
        if hk(2^(-i))<min
            min=hk(2^(-i));
            mini=i;
        end
    end
    
    lambda=2^(-mini);
    if lambda<0.01
        lambda=0.01;
    end
    
    x=x0-lambda*d;
  
  flag=0;
  
  if (norm(x-x0)/norm(x0))<eps
    flag=1;
  end
  
  while flag==0 && (norm(x-x0)/norm(x0))>=eps && k<N
    x0=x;
    F=f(x);
    J=df(x);
    d=J\F;
    gama=1/(cond(J));

  h=@(x) f(x)'*f(x);
  hk=@(tau) h(x-tau*d);
  
    nasao=0;
    j=0;
    while(nasao==0)
        lijevi=hk(2^(-j));
        desni=hk(0)-(2^(-j))*(gama/4)*norm(d)*norm(2*F'*J);
        if lijevi<=desni
            nasao=1;
        end
        j=j+1;
    end
    min=hk(1);
    mini=0;
    for i=0:j
        if hk(2^(-i))<min
            min=hk(2^(-i));
            mini=i;
        end
    end
    
    lambda=2^(-mini);
    if lambda<0.01
        lambda=0.01;
    end
    
    x=x0-lambda*d;
    k=k+1;
    
  end

end