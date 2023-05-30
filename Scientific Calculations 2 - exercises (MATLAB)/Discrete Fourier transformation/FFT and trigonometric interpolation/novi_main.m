f = @(x) exp(-x.^2/4);

n=4;
N=16;

i=0:N-1;
x=(pi*i*2)/N;

bete=FFT(f(x),n)

for i=1:N
  Hf(i)=real(fazni_Horner(x(i),bete));
end

greska=max(abs(Hf-f(x)))

figure()
fplot(f,[0,6],'b-');
hold on;
plot(x,Hf,'ko');
hold off;
legend('egzaktno rjesenje','interpolacijske tocke');

[A,B]=trig_FFT(f(x),n)

h2=(2*pi)/200;
x2=0:h2:2*pi;

for i=1:201
  psievi(i)=gen_horner_trig(x2(i),A,B);
end

figure()
fplot(f,[0,2*pi],'b-');
hold on;
plot(x,Hf,'ko');
plot(x2,psievi,'r-');
hold off;
legend('egzaktno rjesenje','interpolacijske tocke','trigonometrijski polinom');