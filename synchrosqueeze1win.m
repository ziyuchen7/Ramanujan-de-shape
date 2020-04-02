function [tfr, rtfr] = synchrosqueeze1win(tfr, ifd, alpha, h, num_s, fr, HighFreq, fs, ths, varargin)
% Code by Su Li
Smooth = 0;
Reject = 0;
for var_i = 1:length(varargin)
    if strcmp(varargin{var_i}, 'Smooth')
        Smooth = varargin{var_i + 1};
    end
    if strcmp(varargin{var_i}, 'Reject')
        Reject = varargin{var_i + 1};
    end
end

[M, N] = size(tfr);
K = length(-0.5+alpha:alpha:0.5);
if Reject>0
    TH = Reject/fr; %*K/length(h); % should be M

    omega = ifd;
    if num_s ==1
        tfr(abs(omega)>TH/2)=0;
    else
        tfr(abs(omega)>TH/2)=0;
        for kk = 2:num_s
    omega_temp = inf(size(omega));
    omega_size = length(kk+1:kk:M);
    omega_temp(2:omega_size+1,:) = omega(kk+1:kk:M,:);
    tfr(abs(omega_temp)>TH/2)=0;
        end
    end
else
    omega = ifd;
end

tfr = tfr(1:round(HighFreq*fs/fr),:);
omega = round(omega(1:round(HighFreq*fs/fr),:));

[M, N] = size(tfr);
OrigIndex = repmat((1:M)', [1 size(tfr,2)]);
omega(OrigIndex - omega < 1+2*Smooth | OrigIndex - omega > M-2*Smooth)=0;

Ex = mean(sum(abs(tfr)));
Threshold = ths*Ex;	% originally it was 1e-6*Ex
tfr(abs(tfr) < Threshold) = 0;

totLength = size(tfr,1)*size(tfr,2);
new_idx = (1:totLength)'-omega(:);
if Smooth == 0
    rtfr = accumarray([1; new_idx(2:totLength-1); totLength],tfr(:));
else
    SmoothWin = triang(1+2*Smooth)./sum(triang(1+2*Smooth));
    rtfr = accumarray([1; new_idx(2:totLength-1); totLength], SmoothWin(1).*tfr(:));
    for ii = 1:Smooth-1
        rtfr = rtfr + accumarray([1:ii; new_idx(ii+1:totLength-ii)-ii; totLength-ii+1:totLength], SmoothWin(ii+1).*tfr(:));
        rtfr = rtfr + accumarray([1:ii; new_idx(ii+1:totLength-ii)+ii; totLength-ii+1:totLength], SmoothWin(ii+1).*tfr(:));
    end
end
rtfr = [rtfr; zeros(totLength-length(rtfr),1)];

rtfr = reshape(rtfr,size(tfr,1),size(tfr,2));