close all ; clear all ;
% Code to generate figures 1, 3, 4, 7, 8, 9, 12

% basic parameters for STFT
basicTF.win = 1201; %4096;
basicTF.hop = 101; %441;
basicTF.fs = 250;
basicTF.fr = 10/250; % frequency resolution/sampling freq
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here (Hamming and its derivative)
advTF.Smo = 0; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection; 
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 10/250; % highest frequency/sampling freq
advTF.LowFreq = 0.1/250; % lowest frequency/sampling freq
advTF.lpc = 0;
P.num_s = 1; % (harmonic 的根數, 只做簡單探討的話設1就好)
P.num_c = 1;

% parameters for cepstral representation
cepR.g = 0.1; % for generalized cepstrum: 設0 -> log cepstrum; 設2 -> autocorrelation (for single component 可以設大一點 (e.g., 2); for multiple component 經驗值: 0.1~0.3)
cepR.Tc=0; %1E-4; % Global threshold of cepstrum


x=importdata('a17.mat');%44-3
x0 = x(:,5);

Trend = zeros(size(x0)) ;
for ii = 1: length(Trend)
	idx = [max(1,ii-50):min(length(Trend),ii+50)] ;
	Trend(ii) = median(x0(idx)) ;
end

% Detrend step
x = x0 - Trend ; 

x = resample(x,1,4);
x0 = resample(x0,1,4);

t = [1:length(x)]'/basicTF.fs ;

fs = basicTF.fs;
% [b,a] = butter(9,0.01,'high');
% x = filter(b,a,x);
% y=resample(x(:,1),1,10);
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
[tfr0, ~, ~, ~, ~, ~, ~] = CFPH_all(x, basicTF, advTF, cepR, P);

df = (tfrtic(2)-tfrtic(1))*basicTF.fs;
max_freq = 60;
f = max_freq/df + 1;
tfr0 = tfr0(1:f,:);

time_stamp = basicTF.hop/basicTF.fs;
newtfr = abs(tfr0).^0.1;

s = PD_Lasso(newtfr(:,60),100,'Ramanujan',0.01);

% ss1 = PD_Lasso(x,200,'Ramanujan',5,0); % run RPT directly on the signal in time domain
% ss2 = PD_Lasso(x,200,'Ramanujan',20,1);
% ss3 = PD_Lasso(x,200,'Ramanujan',5,2);

stpr1 = RPT_time_period(x,200,'Ramanujan',500,500,101,0); % run short-time RPT directly on the signal in time domain, zeta(p) = 1
stpr2 = RPT_time_period(x,200,'Ramanujan',200,500,101,1); % run short-time RPT directly on the signal in time domain, zeta(p) = p
stpr3 = RPT_time_period(x,200,'Ramanujan',10,500,101,2); % run short-time RPT directly on the signal in time domain, zeta(p) = p^2


if(1)
numper = 100;
np = 3;

tpr1 = tfr_small2large(newtfr,0.01,100); % small to large
tprr1 = tfr0.*tpr1;

tpr2 = tfr_mbest(newtfr,np,numper); % m-best
tprr2 = tfr0.*tpr2;

tpr3 = tfr_bestcor(newtfr,numper,np); % best-correlation
tprr3 = tfr0.*tpr3;

tpr4 = RDS(newtfr,100,'Ramanujan',0.01);
tprr4 = tfr0.*tpr4;

% tpr5 = tfr_L1(newtfr,100,'Ramanujan'); % l^1 minimization
% tprr5 = tfr0.*tpr5;
% 
% tpr6 = tfr_L2(newtfr,100,'Ramanujan'); % l^2 minimization
% tprr6 = tfr0.*tpr6;

[ntpr] = vector_RDS(newtfr,100,'Ramanujan',0.03,1);
ntprr = tfr0.*ntpr;
end


p = 0.999;
figure
plot((1:1000)/250, x0(1:1000),'color','k','LineWidth',2); axis tight; xlabel('time (s)');ylabel('Arbitrary Unit');
set(gca, 'fontsize', 15);axis tight

figure
subplot(3,2,1)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic*basicTF.fs, abs(tfr), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');
%annotation('arrow',[0.2,0.22],[0.56,0.58],'Color','red')
subplot(3,2,2)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic*basicTF.fs, abs(tfr).^0.1, p); axis xy; colormap((1-gray).^(1/4));
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');
subplot(323)
plot((tfrtic(2)-tfrtic(1))*(1:1500)*basicTF.fs, abs(tfr0(:,60)),'color','k');
xlabel('frequency (Hz)'); ylabel('Arbitrary Unit');set(gca, 'fontsize', 15);axis tight
subplot(324)
plot((tfrtic(2)-tfrtic(1))*(1:1500)*basicTF.fs, abs(newtfr(:,60)),'color','k');
xlabel('frequency (Hz)'); ylabel('Arbitrary Unit');set(gca, 'fontsize', 15);axis tight
subplot(3,2,[5,6])
stem(s,'linewidth',3,'color',[0 0 0]);
xlabel('period');axis tight;
set(gca, 'fontsize', 15);

figure
subplot(121)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tceps(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');
subplot(122)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tfrr(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('De-shape STFT');

figure
subplot(121)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tpr4(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Mask');
subplot(122)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tprr4(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('l1 penalized');

figure
subplot(321)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tprr1(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); colorbar%title('small2large');
subplot(322)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tprr2(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); colorbar%title('Mbest');
subplot(323)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tprr3(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); colorbar%title('Best-correlation');
subplot(324)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tprr4(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); colorbar%title('l1 penalized');
subplot(325)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, nptrr(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); colorbar%title('Best-correlation');



figure
subplot(221)
imageSQ(0:time_stamp:time_stamp*148, 1:200, stpr1, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('Period'); %title('STPT 1');
subplot(222)
imageSQ(0:time_stamp:time_stamp*148, 1:200, stpr2, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('Period'); %title('STPT 2');
subplot(223)
imageSQ(0:time_stamp:time_stamp*148, 1:200, stpr3, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('Period'); %title('STPT 3');
