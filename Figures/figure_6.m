clear all; close all;

% Code for figure 6
% Robustness test with respect to gaussian with envelope and additive noise
initstate(3);
N = 301;
x1 = zeros(N,1);
x2 = zeros(N,1);
orig_x1 = zeros(N,1);
orig_x2 = zeros(N,1);


t = 4;
[ch,Dh,tt] = hermfun(2*t+1,1,8);
h = ch(1,:);

amp1 = 3;
amp2 = 3;
T1 = 27;
T2 = 17;
D = 0.2;
for i = T1:T1:N
    orig_x1(i-t:i+t) = amp1.*h';
    d = rand;
    d = 2*D*d-D;
    x1(i-t:i+t) = orig_x1(i-t:i+t).*(1+d);
end

t = 6;
[ch,Dh,tt] = hermfun(2*t+1,1,12);
h = ch(1,:);

for i = T2:T2:N
    orig_x2(i-t:i+t) = amp2.*h';
    d = rand;
    d = 2*D*d-D;
    x2(i-t:i+t) = orig_x2(i-t:i+t).*(1+d);
end

x = x1 + x2;
sigma = 2;
t = (1:N)/N;
env = exp(-t.^2./sigma^2)';

x = x.*env;
orig_x = orig_x1 + orig_x2;
SNR = 5;
ns = randn(length(x),1);%*10^(-1*SNR*norm(x)/20);
ns = ns/norm(ns)*norm(x)*10^(-SNR/20);
x = x + ns;


s = PD_Lasso(x,50,'Ramanujan',0.2,2);
k = 15;
SNR = 20*log10(std(orig_x)/std(ns))
figure
subplot(421)
plot(orig_x1,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(422)
plot(x1,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(423)
plot(orig_x2,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(424)
plot(x2,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(425)
plot(x1+x2,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(426)
plot(env,'color',[0 0 0])
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(427)
plot(x,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(428)
stem(s,'linewidth',3,'color',[0 0 0]);
xlabel('period');xlim([1 50]);set(gca, 'fontsize', k);
