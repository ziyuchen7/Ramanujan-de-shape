clear all; close all;
% Code for figure 5
% Robustness test with respect to envelope
initstate(1);
N = 500;
x1 = zeros(N,1);
x2 = zeros(N,1);
orig_x1 = zeros(N,1);
orig_x2 = zeros(N,1);

amp1 = 5;
amp2 = 8;
T1 = 15;
T2 = 21;
D = 0.3;
for i = T1:T1:N
    orig_x1(i) = amp1;
    d = rand;
    d = 2*D*d-D;
    x1(i) = amp1*(1+d);
end

for i = T2:T2:N
    orig_x2(i) = amp2;
    d = rand;
    d = 2*D*d-D;
    x2(i) = amp2*(1+d);
end
x1(15) = 0;
x2(21) = 0;
x = x1 + x2;
sigma = 1;
t = (1:N)/N;
env = exp(-t.^2./sigma^2)';

x = x.*env;


s = PD_Lasso(x,50,'Ramanujan',0.8,1);
k = 15;

%% small to large
[s1] = small2large(x,0.3,50,1);

%% M-best
[s2] = mbest(x,2,50,1);

%% best correlation
[s3] = bestcorrelation(x,50,2,1);


%%
figure
subplot(621)
stem(orig_x1,'linewidth',3,'color',[0 0 0]);
xlabel('time');axis tight;
set(gca, 'fontsize', k);
subplot(622)
stem(x1,'linewidth',3,'color',[0 0 0]);
xlabel('time');axis tight;
set(gca, 'fontsize', k);
subplot(623)
stem(orig_x2,'linewidth',3,'color',[0 0 0]);
xlabel('time');axis tight;
set(gca, 'fontsize', k);
subplot(624)
stem(x2,'linewidth',3,'color',[0 0 0]);
xlabel('time');axis tight;
set(gca, 'fontsize', k);
subplot(625)
stem(x1+x2,'linewidth',3,'color',[0 0 0]);
xlabel('time');axis tight;
set(gca, 'fontsize', k);
subplot(626)
plot(env,'color',[0 0 0])
xlabel('time');axis tight;
set(gca, 'fontsize', k);
subplot(627)
stem(x,'linewidth',3,'color',[0 0 0]);
xlabel('time');axis tight;
set(gca, 'fontsize', k);
subplot(628)
stem(s,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;
set(gca, 'fontsize', k);
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
