function [s] = bestcorrelation(x,Pmax,n)
% INPUT:
%  x - input signal
%  Pmax - upper bound for period
%  n - number of estimated periods
% OUTPUT:
%  s - strength of period
% Written by Ziyu Chen, Feb 2020

if(size(x,2)>size(x,1))
    x = x';
end

Nmax = Pmax;
A = Create_Dictionary(Nmax,size(x,1),'Ramanujan');

s = 0.*[1:Nmax];

r = x;

for i = 1:n
    energy = 0.*[1:Nmax];
    temp_energy = 0;
    temp_proj = [];
    for j = 1:Nmax
        index = (count_euler(j-1)+1):count_euler(j);
        D = A(:,index);
        proj = D*inv(D'*D)*D' * r;
        energy(j) = proj'*proj;
        if energy(j)>temp_energy
            temp_energy = energy(j);
            temp_proj = proj;
        end
    end
    [a,b] = max(energy);
    s(b) = s(b) + a;
    r = r - temp_proj;
end
    


s(1) = 0;


end

