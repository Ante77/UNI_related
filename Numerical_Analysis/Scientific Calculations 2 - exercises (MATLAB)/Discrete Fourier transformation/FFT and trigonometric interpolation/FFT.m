function b = FFT(f,n)
  
  b=zeros(n,1);
  
  for i=1:2^n
    b(i)=f(rho(i-1,n)+1);
  end
  
  for m=1:n
    for j=1:2^(m-1)
      e=exp((-2*pi*(j-1)*sqrt(-1))/(2^m));
      for q=0:2^m:2^n-1
        u=b(q+j);
        v=b(q+j+2^(m-1))*e;
        b(q+j)=u+v;
        b(q+j+2^(m-1))=u-v;
      end
    end
  end
  
  b=b/(2^n);
  
end