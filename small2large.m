function [x]=small2large(s,thresh,numper)

% Calculate the Periodicity Transform using the
% "Small To Large" algorithm
% Ziyu Chen, Feb 2020
% inputs: s = signal to decompose
%         thresh = P_p must remove at least this percent of 
%                  power in order to be counted
%         numper = max periodicity to look for
% outputs: x = strength of each period
%
% Modified from:
% See Sethares and Staley, "Periodicity Transforms"
% IEEE Trans. Signal Processing, 1999.

A = Create_Dictionary( numper, length(s), 'Ramanujan' );
[nr,nc]=size(s); if nc>nr, s=s'; end
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


end
