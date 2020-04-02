function [tfr0, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH_all(x, basicTF, advTF, cepR, P)
% tfr: STFT; ceps: Cepstrum; tceps: ;tfrr: ;rtfr: ; tfrsq: SST; tfrtic: frequencies in SST
% basic parameters for STFT
% Code by Su Li and Hau-tieng Wu, Ziyu Chen added one more output: tfr0 - the whole STFT
win = basicTF.win;
hop = basicTF.hop;
feat = basicTF.feat;
fs = basicTF.fs;
fr = basicTF.fr;

% advanced parameters for STFT
HighFreq = advTF.HighFreq;
LowFreq = advTF.LowFreq;
num_tap = advTF.num_tap;
win_type = advTF.win_type;
Smo = advTF.Smo;
Rej = advTF.Rej;
ths = advTF.ths;
lpc = advTF.lpc;

num_s = P.num_s;
num_c = P.num_c;

h = tftb_window(win, win_type);
Dh = dwindow(h);
Dho = dwindow(Dh);

h = [h Dh];
Dh = [Dh Dho];

% parameters for cepstral representation
% med_num = cepR.med_num; %11;
g = cepR.g; %0.06;
Tc = cepR.Tc;

MT = num_tap;
for ii = 1:MT
fprintf(['ConceFT total: ',num2str(MT),'; now: %4d\n'], ii) ;

if MT == 1
    rv = [1 0];
else
    rv = randn(1, 2) ;  %rv(2)=0;
    rv = rv ./ norm(rv) ;
end
rh = rv * h';
rDh = rv * Dh';
h2 = rh'; Dh2 = rDh';

% x = x(1000:5000);
% a = lpc(x,150);
% est_x = filter(a,1,x);

[tfr, ifd, tfrtic, Stime] = STFT_IFD_fast(x, fr/fs, hop, h2, Dh2);

tfr = abs(tfr); % STFT
% tfr = [diff(tfr,[],2) zeros(size(tfr,1),1)];
% tfr(tfr<0)=0;

%%% ceps
if g~=0
    ceps = 2.*real(ifft(abs(tfr).^g,2*size(tfr,1),1));
    bnd = floor(1/HighFreq);
    ceps(bnd+1:end-bnd,:)=0;
    flo = (real(fft(ceps))).^(1/g); flo = flo(1:size(tfr,1),:);
else
    ceps = 2.*real(ifft(log(abs(tfr)),2*size(tfr,1),1));
    bnd = floor(1/HighFreq);
    ceps(bnd+1:end-bnd,:)=0;
    flo = (real(fft(ceps))).^(1/g); flo = flo(1:size(tfr,1),:);
end

if lpc>0
    [tfr_lpc, ~, ~, ~] = STFT_IFD_lpc_fast(x, fr/fs, hop, h2, Dh2, lpc);
    [ceps, tceps] = cepstrum_convert(tfr_lpc, tfrtic, g, fs, Tc, num_c, HighFreq, LowFreq);
else
    [ceps, tceps] = cepstrum_convert(tfr, tfrtic, g, fs, Tc, num_c, HighFreq, LowFreq);
end
% tceps(tceps>0)=1;
tfr0 = tfr; 
tfrr = tfr0.*tceps;
tfrr(tfrr<0)=0;
tfrr(1:round(LowFreq*fs/fr),:)=0;

tfrr = tfrr(1:round(HighFreq*fs*num_s/fr),:);
ifd = ifd(1:round(HighFreq*fs*num_s/fr),:);

[tfr2, tfr3] = synchrosqueeze1win(tfrr, ifd, fr/fs, h2, num_s, fr, HighFreq, fs, ths, 'Smooth', Smo, 'Reject', Rej);
[~, tfrsq] = synchrosqueeze1win(tfr0, ifd, fr/fs, h2, num_s, fr, HighFreq, fs, ths, 'Smooth', Smo, 'Reject', Rej);

if MT == 1
    if strcmp(feat,'STFT')
        rtfr = tfr2;
    elseif strcmp(feat,'SST11')
        rtfr = tfr3;
    end
else % MT>1, no STFT
    if ii == 1
        if strcmp(feat,'STFT')
            rtfr = tfr2;
        elseif strcmp(feat,'SST11')
            rtfr = tfr3;
        end
    else
        if strcmp(feat,'STFT')
            rtfr = rtfr+tfr2;
        elseif strcmp(feat,'SST11')
            rtfr = rtfr + tfr3;
        end
    end
end
end

tfrr = tfrr(1:round(HighFreq*fs/fr),:);
tceps = tceps(1:round(HighFreq*fs/fr),:);
tfr = tfr(1:round(HighFreq*fs/fr),:);

start_idx = floor(Stime/fs/(hop/fs));
tfr = [zeros(size(rtfr,1), start_idx) tfr];
tfr0 = [zeros(size(tfr0,1), start_idx) tfr0];
tfrr = [zeros(size(rtfr,1), start_idx) tfrr];
tceps = [zeros(size(rtfr,1), start_idx) tceps];
rtfr = [zeros(size(rtfr,1), start_idx) rtfr];
tfrtic = tfrtic(1:size(rtfr,1));
% rtfr = medfilt2(rtfr,[1 med_num]);

