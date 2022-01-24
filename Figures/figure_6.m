clear all; close all;
% Robustness test with respect to random translations of peaks
initstate(1);
N = 300;
x = zeros(N,1);
y = zeros(N,1);

T = 19;
D = 1;
amp = 7;

for i = T:T:N
    j = i;
    d = rand;
    d = round(2*D*d-D);
    j = j + d;
    x(j) = amp;
end

for i = T:T:N
    y(i) = amp;
end

% T2 = 7;
% for i = T2:T2:N
%     j = i;
%     d = rand;
%     d = round(2*D*d-D);
%     j = j + d;
%     x(j) = x(j) + 7;
% end
% 
% for i = T2:T2:N
%     y(i) = y(i) + 10;
% end

%s = Regularized_Strength_vs_Period_L1_figure(x,100,'Ramanujan',0.001);
%s = Normalized_Regularized_Strength_vs_Period_L1_figure(x,100,'Ramanujan',0.1);
s = PD_Lasso(x,50,'Ramanujan',1,2);

%% small to large
[s1] = small2large(x,0.5,50,2);
%% M-best
[s2] = mbest(x,1,50,2);
%% best correlation
[s3] = bestcorrelation(x,50,1,2);

%%
k = 15;
figure()
subplot(421)
plot(y,'color','k');axis tight;
xlabel('time');set(gca, 'fontsize', k);
subplot(422)
plot(x,'color','k');axis tight;
xlabel('time');set(gca, 'fontsize', k);
subplot(423)
plot(x,'color','r');axis tight;
hold on 
plot(y,'color','k');axis tight;
xlabel('time');
set(gca, 'fontsize', k);
subplot(424)
stem(s,'linewidth',3,'color',[0 0 0]);
xlabel('Period');axis tight;
set(gca, 'fontsize', k);
subplot(425)
stem(s1,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;
set(gca, 'fontsize', k);
subplot(426)
stem(s2,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;
set(gca, 'fontsize', k);
subplot(4,2,7.5)
stem(s3,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;
set(gca, 'fontsize', k);
