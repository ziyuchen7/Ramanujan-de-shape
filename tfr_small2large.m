function [tpr]=tfr_small2large(tfr,thresh,numper)

% Time-frequency representation by Small-to-Large periodicity transform algorithm
% Written by Ziyu Chen, Feb 2020
% syntax: [x]=small2large(s,thresh,numper)
% inputs: tfr = input time frequency representation
%         thresh = P_p must remove at least this percent of 
%                  power in order to be counted
%         numper = max periodicity(frequency)
% outputs: tpr = result time frequency representation
%
% Modified from:
% See Sethares and Staley, "Periodicity Transforms"
% IEEE Trans. Signal Processing, 1999.

tpr = zeros(size(tfr,1),size(tfr,2));
A = Create_Dictionary( numper, size(tfr,1), 'Ramanujan' );
for l = 2:size(tfr,2)
    s = tfr(:,l);
    x = zeros(numper,1);
    snorms=periodnorm(s);
    if nargin==2
        numper=length(s)/2;
    end
    
    for period=1:numper
        D = A(:,count_euler(period-1)+1:count_euler(period));
        bas = D*inv(D'*D)*D'*s;
        r=s-bas;
        normimp = periodnorm(bas)/snorms;
        if normimp>thresh
            s=r;
            x(period) = periodnorm(bas)^2;
        end
    end
    x(1) = 0;
    tpr(2:numper+1,l) = x;
end

end
