function [x,y,s,k]=odj_gadjanje_2rru_en(f,dfy,dfdy,a,b,alfa,beta,s0,n) 
  
  x=zeros(1,n+1);
  y=zeros(4,n+1);
  y0=[alfa; s0; 0; 1];
  epsilon=1e-14;
  s=s0;
  k=0;
  
  while abs(y(1,n+1)-beta)>=epsilon || k==0
    y0=[alfa; s; 0; 1];
    
    F=@(x,y) [f(x,y); dfy(x,y); dfdy(x,y)];    
    [x,y]=odj_rk4v(F,a,b,y0,n);
    
    s=s-(y(1,n+1)-beta)/y(3,n+1);
    k=k+1;
  end
  k=k-1;
  
end