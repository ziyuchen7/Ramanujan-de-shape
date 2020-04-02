function [ceps0, tceps] = cepstrum_convert(tfr, tfrtic, g, fs, Tc, num_s, HighFreq, LowFreq)
% Code written by Su Li for inverse STCT 
if g~=0
    ceps = real(ifft(abs(tfr).^g,2*size(tfr,1),1));
else
    ceps = real(ifft(log(abs(tfr)),2*size(tfr,1),1));
end
for mi=1:size(ceps,2)
    tra = min(find(ceps(3:round(1/LowFreq),mi)<0))-1+3;
    if tra>=1
        ceps(1:max([round(1/HighFreq) tra]),mi)=0;
    else
    end
end

ceps(1:round(1/HighFreq),:)=0; 
% ceps(ceps<Tc)=0; %ceps(ceps>0)=1;
% % 
% if num_s>1
% for kk = 2:num_s
%     ceps_temp = zeros(size(ceps));
%     ceps_size = length(kk+1:kk:size(ceps,1));
%     ceps_temp(2:ceps_size+1,:) = ceps(kk+1:kk:size(ceps,1),:);
%     ceps = ceps.*ceps_temp;
% end
% end

UpSample = 20 ;

ceps = ceps(1:round(1/LowFreq),:);
ceps0 = ceps;
ceps =  interp1(1:size(ceps,1), ceps, 1:1/UpSample:size(ceps,1));
tceps = zeros(length(tfrtic), size(ceps,2));
	% cepstrum quefrency scale
freq_scale = UpSample.*fs./(1:size(ceps,1)-1);
empt = [];


for ii = 2:length(tfrtic)-1
		% index in quefency
    p_index = find(freq_scale > (tfrtic(ii-1)+tfrtic(ii))*fs/2 & freq_scale <= (tfrtic(ii+1)+tfrtic(ii))*fs/2);
    if isempty(p_index)
        empt = [empt; ii];
    else
        %tceps(ii,:)=sum(ceps(p_index,:),1);
		weight = abs(1./p_index) ;
        tceps(ii,:)=sum(diag(weight)*ceps(p_index,:),1);
    end
end
tceps(tceps<Tc)=0;
tceps_new = tceps;

	% how many peaks in cepstrum (num_s)
if num_s>1
    for fi = round(LowFreq*fs*num_s)+1:size(tceps,1)
        for kk = 2:num_s
            tceps_new(fi,:)=tceps(fi,:).*tceps(max([1 round(fi/kk)]),:);
        end
    end
% for kk = 2:num_s
%     tceps_temp = zeros(size(tceps));
%     tceps_size = length(kk+1:kk:size(tceps,1));
%     tceps_temp(2:ceps_size+1,:) = ceps(kk+1:kk:size(ceps,1),:);
%     ceps = ceps+ceps_temp;
% end
end
tceps = tceps_new;
% ceps_low = interp1(2:50, ceps(2:50,:), 2:0.001:50);
% freq_scale_low=1000.*fs./(1:49);
% for ii=2:length(tfrtic)-1
%     p_index = find(freq_scale_low > (tfrtic(ii-1)+tfrtic(ii))*fs/2 & freq_scale_low < (tfrtic(ii+1)+tfrtic(ii))*fs/2);
%     if isempty(p_index)
% %         empt = [empt; ii];
%     else
%         tceps(ii,:)=sum(ceps_low(p_index,:),1);%./(ii);
%     end
% end
% 
% ceps_mid = interp1(51:200, ceps(51:200,:), 51:0.01:200);
% freq_scale_mid = 100.*fs./(50:199);
% for ii=2:length(tfrtic)-1
%     p_index = find(freq_scale_mid > (tfrtic(ii-1)+tfrtic(ii))*fs/2 & freq_scale_mid < (tfrtic(ii+1)+tfrtic(ii))*fs/2);
%     if isempty(p_index)
% %         empt = [empt; ii];
%     else
%         tceps(ii,:)=sum(ceps_mid(p_index,:),1);%./(ii);
%     end
% end
% 
% 
% ceps_high = interp1(201:size(ceps,1), ceps(201:end,:), 201:0.1:size(ceps,1));
% freq_scale_high = 10.*fs./(200:size(ceps,1)-1);
% for ii=2:length(tfrtic)-1
%     p_index = find(freq_scale_high > (tfrtic(ii-1)+tfrtic(ii))*fs/2 & freq_scale_high < (tfrtic(ii+1)+tfrtic(ii))*fs/2);
%     if isempty(p_index)
% %         empt = [empt; ii];
%     else
%         tceps(ii,:)=sum(ceps_high(p_index,:),1);%./(ii);
%     end
% end

%tceps(tceps>0)=1;
