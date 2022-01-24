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

%% small to large
[s1] = small2large(x,0.3,50,2);
%% M-best
[s2] = mbest(x,2,50,2);
%% best correlation
[s3] = bestcorrelation(x,50,2,2);

%%
k = 15;
SNR = 20*log10(std(orig_x)/std(ns))
figure
subplot(621)
plot(orig_x1,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(622)
plot(x1,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(623)
plot(orig_x2,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(624)
plot(x2,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(625)
plot(x1+x2,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(626)
plot(env,'color',[0 0 0])
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(627)
plot(x,'color',[0 0 0]);
xlabel('time');axis tight;set(gca, 'fontsize', k);
subplot(628)
stem(s,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;set(gca, 'fontsize', k);
subplot(629)
stem(s1,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;
set(gca, 'fontsize', k);
subplot(6,2,10)
stem(s2,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;
set(gca, 'fontsize', k);
subplot(6,2,11.5)
stem(s3,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;
set(gca, 'fontsize', k);