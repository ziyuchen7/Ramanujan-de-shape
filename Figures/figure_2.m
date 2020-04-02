close all; clear all;


initstate(1) ;
dd = zeros(1e5+1,1) ;
dd(1:1000:end) = 1 ;
aa = exp(-([-2e4:8e4]/1e5).^2) ;
[h,dh] = hermf(201,1,6) ;
h = h ./ max(h) ;

Hz = 1000 ;
t = [1:100001]' / Hz ;
x1 = 3*conv(aa'.*dd,h,'same') ; % 1Hz
x1(60001:end) = 0 ;

ff = abs(cumsum(randn(size(x1)))) ; IF2 = ff./(max(abs(ff))/2) + pi/2 ;
IF2 = smooth(IF2, 10000) ;
phi = cumsum(IF2) ./ Hz ;
AM2 = smooth(abs(cumsum(randn(size(x1)))./Hz) + 1, 20000) ;
AM2 = AM2 ./ max(AM2) + .9 ;
gg = mod(phi,1);
[a,b] = findpeaks(gg);
b = [1; b; 2*b(end)-b(end-1)] ;
s2 = zeros(size(phi)) ;
for ii = 1: length(b)-1
    idx = b(ii):b(ii+1) ;
    s2(idx) = (idx-b(ii)) ./ (b(ii+1)-b(ii)+1) ;
end
x2 = 5 .* s2(1:length(AM2)) ;
s = @(x) cos(2*pi*x) + 2*cos(4*pi*x) + 3*cos(6*pi*x) + 4*cos(8*pi*x) + ...
    7*cos(10*pi*x) + 10*cos(12*pi*x) + 8*cos(14*pi*x) + 6*cos(16*pi*x) + ...
    3*cos(18*pi*x) + cos(20*pi*x);
x3 = s(phi);
x = x3 ;
x = x - mean(x);
noise = randn(size(x)) * 0.5 ;
snrdb = 20*log10(std(x)/std(noise))
x = x(1:20:end); %+ noise(1:20:end) ;
t = t(1:20:end) ;
Hz = Hz / 20 ;

basicTF.win = 300; %4096;
basicTF.hop = 21; %441;
basicTF.fs = Hz;
basicTF.fr = 0.02;
basicTF.feat = 'SST11';
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = 'Gauss'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
advTF.Smo = 0;
advTF.Rej = 0;
advTF.ths = 1E-9;
advTF.HighFreq = 16/50;
advTF.LowFreq = 0.1/50;
advTF.lpc = 0;
cepR.g = 0.1; % g(cc);%0.06;
cepR.Tc=0;
P.num_s = 1;
P.num_c = 1;

[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
% rtfr = medfilt2(rtfr,[1 3]);
time_stamp = basicTF.hop/basicTF.fs;


if(1)
p = 0.999;
figure
subplot(2,1,1)
plot((1:1000)/250, x(1:1000),'color','k'); axis tight; xlabel('time (s)');ylabel('Arbitrary Unit');
set(gca, 'fontsize', 15);axis tight
subplot(2,1,2)
imageSQ(0:time_stamp:time_stamp*(size(tfr,2)-1), tfrtic*basicTF.fs, tfr, p); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');
end