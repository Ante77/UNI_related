function [f] = fja(x,y,a)
    
    f=x.^2 +a(1)*x.*y+a(2)*y.^2+a(3)*x+a(4)*y+a(5);
    
end