function [min,it]=Powell_algorithm2(f,x0,epsilon)
  dim=max(size(x0));
  xn=x0;
  S=zeros(dim,dim+1);
  S(1:dim,1:dim)=eye(dim);
  alfe=zeros(dim,1);
  it=1;
  X(:,1)=x0;
  xpom=zeros(dim,dim+1);
  
  while norm(xn-x0)>=epsilon || it==1
    x0=xn;
    xpom(:,1)=x0;
    for i=1:dim
      fja=@(alfa) f(xpom(:,i)+alfa*S(:,i));
      alfe(i)=fminbnd(fja,-100,100);
      xpom(:,i+1)=xpom(:,i)+alfe(i)*S(:,i);
      S(:,i)=S(:,i+1);
    end
    spravi=xpom(:,dim+1)-xpom(:,1);
    fjaprava=@(alfa) f(xpom(:,dim+1)+alfa*spravi);
    alfaminpravi=fminbnd(fjaprava,-100,100);
    S(:,dim)=spravi;
    Sovi=S(:,1:dim)
    xpom
    alfe
    alfaminpravi
    xn=xpom(:,dim+1)+alfaminpravi*spravi
    it=it+1;
  end
  min=xn;
  it=it-1;
end
