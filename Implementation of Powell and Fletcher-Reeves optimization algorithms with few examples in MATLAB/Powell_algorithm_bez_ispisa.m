function [min,it]=Powell_algorithm_bez_ispisa(f,x0,epsilon)
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
      alfe(i)=Quadratic_search_algorithm(-1000,1000,fja,epsilon);
      xpom(:,i+1)=xpom(:,i)+alfe(i)*S(:,i);
      S(:,i)=S(:,i+1);
    end
    spravi=xpom(:,dim+1)-xpom(:,1);
    fjaprava=@(alfa) f(xpom(:,dim+1)+alfa*spravi);
    alfaminpravi=Quadratic_search_algorithm(-1000,1000,fjaprava,epsilon);
    S(:,dim)=spravi;
    xn=xpom(:,dim+1)+alfaminpravi*spravi;
    it=it+1;
  end
  min=xn;
  it=it-1;
end