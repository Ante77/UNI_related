function [xmin]=Quadratic_search_algorithm(x0,x2,f,epsilon)
  x1=(x0+x2)/2;
  h=x1-x0;
  while h>=epsilon
    if f(x0)<f(x1)
      x2=x1;
      x1=x0;
      x0=x1-h;
    end
    if f(x2)<f(x1)
      x0=x1;
      x1=x2;
      x2=x1+h;
    end
    xN=x1+(h*(f(x0)-f(x2)))/(2*(f(x0)-2*f(x1)+f(x2)));
    if f(xN)<f(x1)
      x1=xN;
    end
    h=h/2;
    x0=x1-h;
    x2=x1+h;
  end
  xmin=x1;
end