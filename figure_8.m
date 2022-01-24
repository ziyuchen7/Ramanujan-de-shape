clear all; close all;
initstate(1)
Hz = 200;
t = 0:1/Hz:60;
N = length(t);

am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am1 = 1 + am1 ./ max(abs(am1))/10 ;
am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am2 = 1.2 + am2 ./ max(abs(am2))/10 ;


    %% the instantaneous IF of the simulated signal
if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
if1 = 2 + if1 ./ max(abs(if1))/2 ;%2
if2 = smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess') ;
if2 = 3.5 + if2 ./ max(abs(if2))/2 ;%4
 
    %% the instantaneous phase of the simulated signal
phi1 = cumsum(if1) / Hz ; 
phi2 = cumsum(if2) / Hz ; 


N1 = 21;
signal1 = zeros(N1,length(t));
for i = 1:N1
    if mod(i,2)==1
        a = 1;
    else a = 0.3;   
    end
    signal1(i,:) = a * cos(2*pi*i*phi1);
end
x1 = am1'.*sum(signal1,1);

N2 = 20;
signal2 = zeros(N2,length(t));
for i = 1:N2
    if mod(i,2)==1
        a = 1.2;
    else 
        a = 0.3;   
    end
    signal2(i,:) = a * cos(2*pi*i*phi2);
end
x2 = am2'.*sum(signal2,1);

x = x1 + x2;

basicTF.win = 684; %1096;
basicTF.hop = 40; %20;
basicTF.fs = Hz;
basicTF.fr = 0.05;
basicTF.feat = 'SST11';
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = 'Gauss'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
advTF.Smo = 0;
advTF.Rej = 0;
advTF.ths = 1E-9;
advTF.HighFreq = 20/Hz;
advTF.LowFreq = 0.1/Hz;
advTF.lpc = 0;
cepR.g = 0.1; % g(cc);%0.06;
cepR.Tc=0;
P.num_s = 1;
P.num_c = 1;

[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
[tfr0, ~, ~, ~, ~, ~, ~] = CFPH_all(x, basicTF, advTF, cepR, P);

df = (tfrtic(2)-tfrtic(1))*basicTF.fs;
max_freq = Hz/2;
f = max_freq/df + 1;
tfr0 = tfr0(1:f,:);
time_stamp = basicTF.hop/basicTF.fs;
newtfr = abs(tfr0).^cepR.g;

%%
tpr = RDS(newtfr,100,'Ramanujan',0.03,2);
tprr = abs(tfr0).*tpr;

[ntpr] = vector_RDS(newtfr,100,'Ramanujan',0.09,1,2);
ntprr = abs(tfr0).*ntpr;

p = 0.99;
figure()
subplot(4,2,1)
plot((2000:10000)/Hz, x1(2000:10000),'color','k','LineWidth',2); xlabel('time (s)');ylabel('Arbitrary Unit');xlim([20 24]);
set(gca, 'fontsize', 18);

subplot(4,2,2)
plot((2000:10000)/Hz, x2(2000:10000),'color','k','LineWidth',2); xlabel('time (s)');ylabel('Arbitrary Unit');xlim([20 24]);
set(gca, 'fontsize', 18);

subplot(4,2,3)
plot((2000:10000)/Hz, x(2000:10000),'color','k','LineWidth',2); xlabel('time (s)');ylabel('Arbitrary Unit');xlim([20 24]);
set(gca, 'fontsize', 18);

subplot(424)
imageSQ((1:size(tfr,2))*time_stamp, tfrtic*basicTF.fs, tfr(:,:), p); axis xy; colormap(1-gray);xlim([10 50]); 
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Short-time Fourier transform');

subplot(425)
imageSQ((1:size(tfr,2))*time_stamp, tfrtic(1:101)*basicTF.fs, tfrr(1:101,:), p); axis xy; colormap(1-gray); xlim([10 50]); 
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('de-shape STFT');

subplot(426)
imageSQ((1:size(tfr,2))*time_stamp, tfrtic(1:101)*basicTF.fs, tprr(1:101,:), p); axis xy; colormap(1-gray); xlim([10 50]); 
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('RDS');

subplot(4,2,7.5)
imageSQ((1:size(tfr,2))*time_stamp, tfrtic(1:101)*basicTF.fs, ntprr(1:101,:), p); axis xy; colormap(1-gray); xlim([10 50]); 
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('vRDS');
