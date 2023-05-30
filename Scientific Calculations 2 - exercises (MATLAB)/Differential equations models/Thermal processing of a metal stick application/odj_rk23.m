function [x,y] = odj_rk23 (f, a, b, y0, tol)

x(1)=a;
y(1)=y0;
Cy(1)=y0;

omega1=214/891;
omega2=1/33;
omega3=650/891;
Comega1=533/2106;
Comega2=0;
Comega3=800/1053;
Comega4=-1/78;

h(1)=b-a;
i=1;

while x(i)<b
  
    if i>1
      h(i)=0.9*h(i-1)*(tol/(abs(y(i)-Cy(i))))^(1/3);
    end
    if x(i)+h(i)>b
      h(i)=b-x(i);
    end
  
    k1=f(x(i),y(i));
    k2=f(x(i)+1/4*h(i),y(i)+h(i)*1/4*k1);
    k3=f(x(i)+27/40*h(i),y(i)+h(i)*((-189/800)*k1+729/800*k2));
    k4=f(x(i)+h(i),y(i)+h(i)*(214/891*k1+1/33*k2+650/891*k3));

    y(i+1)=y(i)+h(i)*(omega1*k1+omega2*k2+omega3*k3);
    Cy(i+1)=y(i)+h(i)*(Comega1*k1+Comega2*k2+Comega3*k3+Comega4*k4);
    
    while abs(y(i+1)-Cy(i+1))>tol
      h(i)=0.9*h(i)*(tol/(abs(y(i+1)-Cy(i+1))))^(1/3);
      
      if x(i)+h(i)>b
        h(i)=b-x(i);
      end
      
      k1=f(x(i),y(i));
      k2=f(x(i)+1/4*h(i),y(i)+h(i)*1/4*k1);
      k3=f(x(i)+27/40*h(i),y(i)+h(i)*((-189/800)*k1+729/800*k2));
      k4=f(x(i)+h(i),y(i)+h(i)*(214/891*k1+1/33*k2+650/891*k3));

      y(i+1)=y(i)+h(i)*(omega1*k1+omega2*k2+omega3*k3);
      Cy(i+1)=y(i)+h(i)*(Comega1*k1+Comega2*k2+Comega3*k3+Comega4*k4);
      end
    
    x(i+1)=x(i)+h(i);
    i=i+1;

      end

    end