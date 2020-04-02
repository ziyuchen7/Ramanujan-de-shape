close all ; clear all ;

% Code to generate figure 10 and 11
% basic parameters for STFT
basicTF.win = 1201; %4096;
basicTF.hop = 51; %441;
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
P.num_s = 1; 
P.num_c = 1;

% parameters for cepstral representation
cepR.g = 0.1; % for generalized cepstrum
cepR.Tc=0; %1E-4; % Global threshold of cepstrum

x=importdata('a31.mat');
x0 = x(:,3);
[a,b] = findpeaks(x0,'MinPeakHeight',0.4) ;
[irr] = interp1((b(2:end)+b(1:end-1))/2, 1000./(b(2:end)-b(1:end-1)), 1:length(x), 'cubic') ;

Trend = zeros(size(x0)) ;
for ii = 1: length(Trend)
	idx = [max(1,ii-50):min(length(Trend),ii+50)] ;
	Trend(ii) = median(x0(idx)) ;
end

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
max_freq = 49;
f = max_freq/df + 1;
tfr0 = tfr0(1:f,:);
%tfr0 = [flip(tfr0(2:end,:),1);tfr0];
time_stamp = basicTF.hop/basicTF.fs;


if(1)
newtfr = abs(tfr0).^0.1;
[tpr1] = RDS(newtfr,100,'Ramanujan',0.01,2);
[tpr2] = RDS(newtfr,100,'Ramanujan',0.03,2);
tprr1 = tfr0.*tpr1;
tprr2 = tfr0.*tpr2;

[ntpr1] = vector_RDS(newtfr,100,'Ramanujan',0.03,1,2);
ntprr1 = tfr0.*ntpr1;
[ntpr2] = vector_RDS(newtfr,100,'Ramanujan',0.09,1,2);
ntprr2 = tfr0.*ntpr2;
end

if(1)
p = 0.999;
figure
subplot(5,2,[1,2])
plot((1:2500)/250, x0(1:2500), 'color','k'); axis tight; xlabel('time (s)');ylabel('Arbitrary Unit');
set(gca, 'fontsize', 15);axis tight
subplot(5,2,[3,4])
plot((1:2500)/250, x(1:2500), 'color','k'); axis tight; xlabel('time (s)');ylabel('Arbitrary Unit');
set(gca, 'fontsize', 15);axis tight
subplot(5,2,5)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic*basicTF.fs, tfr, p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');
subplot(526)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tfrr(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Deshape STFT');
subplot(527)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tprr1(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Period method');
subplot(528)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, tprr2(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Period method');
subplot(5,2,9)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, ntprr1(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Period method');
subplot(5,2,10)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic(1:101)*basicTF.fs, ntprr2(1:101,:), p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Period method');
end