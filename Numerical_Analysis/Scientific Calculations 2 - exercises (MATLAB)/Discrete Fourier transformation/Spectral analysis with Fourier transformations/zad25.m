clear
clc

g=@(t) 0.2*cos(2*pi*20.*t)+0.35*sin(2*pi*50.*t)+0.3*sin(2*pi*70.*t);

N=512;
n=9;
T=2;
dt=T/N;
dfi=1/T;

t=0:dt:(T-dt);
gk=g(t);
gkp=gk+0.5*randn(size(t));

figure()
subplot(2,1,1);
plot(gk,'b-');
xlabel('t');
ylabel('y');
title('Uzorak signala bez suma');
subplot(2,1,2);
plot(gkp,'b-');
xlabel('t');
ylabel('y');
title('Uzorak signala sa sumom');

bete=FFT(gk,n);
betep=FFT(gkp,n);

for i=1:max(size(bete))
    dom_bete(i)=i/T;
end

figure()
subplot(2,1,1);
stem(dom_bete,abs(bete),'r-o');
xlabel('\Phi_k');
ylabel('|\beta_k|');
title('|\beta_k| bez suma');
subplot(2,1,2);
stem(dom_bete,abs(betep),'r-o');
xlabel('\Phi_k');
ylabel('|\beta_k|');
title('|\beta_k| sa sumom');

[a,b]=trig_FFT(gk,n);
[ap,bp]=trig_FFT(gkp,n);

for i=1:max(size(a))
    dom_a(i)=i/T;
end
for i=1:max(size(b))
    dom_b(i)=i/T;
end

figure()
subplot(2,1,1);
stem(dom_a,a,'r-o');
xlabel('\Phi_k');
ylabel('A_k');
title('A_k bez suma');
subplot(2,1,2);
stem(dom_b,b,'r-o');
xlabel('\Phi_k');
ylabel('B_k');
title('B_k bez suma');

figure()
subplot(2,1,1);
stem(dom_a,ap,'r-o');
xlabel('\Phi_k');
ylabel('A_k');
title('A_k sa sumom');
subplot(2,1,2);
stem(dom_b,bp,'r-o');
xlabel('\Phi_k');
ylabel('B_k');
title('B_k sa sumom');