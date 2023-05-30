function [x,y,s,k]=odj_gadjanje_2rru(f,a,b,alfa,beta,s0,n)
  
  x=zeros(1,n+1);
  y=zeros(3,n+1);
  y0=[s0;0;s0];
  epsilon=1e-14;
  s=s0;
  k=0;
  
  while abs(y(1,n+1)-beta)>=epsilon || k==0
    y0=[s;0;s];
    [x,y]=odj_rk4v(f,a,b,y0,n);
    
    deltas=sqrt(eps)*s;
    F=y(1,n+1)-beta;
    y0=[s+deltas;0;s+deltas];
    [x,y2]=odj_rk4v(f,a,b,y0,n);
    DF=(y2(1,n+1)-y(1,n+1))/deltas;
    
    s=s-F/DF
    k=k+1;
  end
  k=k-1;
  
end