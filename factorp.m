function allfactors=factorp(n)

% Find all factors of an integer n
%
% syntax: allfactors=factor(n)
% input:  n = integer to factor
% output: allfactors = vector containing all factors
%

allfactors=[];
for i=1:n/2
  if n/i==round(n/i)
    allfactors=[allfactors, i];
  end
end
