function [x]=mbest(s,m,n)

% Calculate the Periodicity Transform using the
% "M-Best" Algorithm
% Ziyu Chen, Feb 2020

% inputs:  s = signal to decompose
%          m = number of desired periods
%          n = longest periodicity to look for (optional: default=length(s)/3)
% outputs: x = strength of each period
%
% Modified from:
% See Sethares and Staley, "Periodicity Transform"
% IEEE Trans. Signal Processing, 1999.

[nr,nc]=size(s);
if nc>nr, s=s'; end
x = zeros(n,1);
snorms=periodnorm(s); 
if nargin==1, m=10; n=length(s)/3; end
if nargin==2, n=length(s)/3; end
listp=zeros([1,m]); listnorm=zeros([1,m]); listbas=zeros([length(s),m]);
if m>1, listchange=1;else listchange=0;end

A = Create_Dictionary(n,length(s),'Ramanujan');
% build initial list of m best periodicities

for i=1:m
    maxnorm=0; maxp=0; maxbas=zeros([length(s),1]);
    for p=1:n
        D = A(:,count_euler(p-1)+1:count_euler(p));
        bas = D*inv(D'*D)*D'*s;
        if periodnorm(bas)>maxnorm
            maxp=p; maxnorm=periodnorm(bas); maxbas=bas;
        end
    end
    listp(i)=maxp; listnorm(i)=maxnorm; listbas(:,i)=maxbas;
    s=s-maxbas;
end

% step 2: decompose basis elements and residuals

while listchange==1
  i=1;
  while i<m
    listchange=0;
    maxnorm=0;
	fact=factorp(listp(i));
    for nf=1:length(fact)
	  p=fact(nf);
      B = A(:,count_euler(p-1)+1:count_euler(p));
      bas= B*inv(B'*B)*B'*listbas(:,i);
      if periodnorm(bas)>=maxnorm
	    maxp=p; maxnorm=periodnorm(bas); maxbas=bas;
      end
    end
    E = A(:,count_euler(maxp-1)+1:count_euler(maxp));
	xbigq=E*inv(E'*E)*E'*listbas(:,i); xsmallq=listbas(:,i)-xbigq;
	nbigq=periodnorm(xbigq); nsmallq=periodnorm(xsmallq);
	minq=min(listnorm); ptwice=length(find(listp==maxp));
	if nsmallq+nbigq>listnorm(m)+listnorm(i) & nsmallq>minq & nbigq>minq & ptwice==0
	  listchange=1;
	  listnorm=[listnorm(1:i-1),nbigq,nsmallq,listnorm(i+1:m-1)];
	  listp=[listp(1:i-1),maxp,listp(i:m-1)];
	  if i>1
	    listbas=[listbas(:,1:i-1);maxbas;xsmallq;listbas(:,i+1:m-1)];
	  else
	    listbas=[maxbas xsmallq listbas(:,2:m-1)];
	  end
	else
	  i=i+1;
	end
  end
end

periods=listp;
powers=listnorm/snorms;
basis=listbas;
x(periods) = powers;

end
