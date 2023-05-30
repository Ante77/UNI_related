function [x,k,rezici] = sor(A,b,x0,tol,omega)
  normb=norm(b);
  x=x0;
  n=max(size(x));
  k=1;
  
  rezici(1)=norm(A*x0-b)/normb;
  while(rezici(k)>tol)
    k=k+1;
    for i=1:n
        x(i)=x(i)+omega*(b(i)-A(i,1:n)*x(1:n))/A(i,i);
    end
    rezici(k)=norm(A*x-b)/normb;
  end
  k=k-1;

end