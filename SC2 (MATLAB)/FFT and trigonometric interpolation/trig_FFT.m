function [A,B] = trig_FFT(f,n)
  
  bete=FFT(f,n);
  N=2^n;
  M=N/2;
  
  A=zeros(M+1,1);
  B=zeros(M,1);
  
  A(1)=2*bete(1);  
  for h=2:M
    A(h)=bete(h)+bete(N-h+2);
    B(h)=(bete(h)-bete(N-h+2))*sqrt(-1);
  end
  
  A(M+1)=2*bete(M+1);
  
end