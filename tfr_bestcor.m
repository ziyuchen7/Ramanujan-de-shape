function [tpr] = tfr_bestcor(tfr,Pmax,n)
% Time-frequency representation by Best-correlation periodicity transform algorithm
% INPUT:
%  tfr - input TPR
%  Pmax - upper bound for period(frequency)
%  n - number of estimated periods(frequencies)
% OUTPUT:
%  tpr - result time frequency representation
% Written by Ziyu Chen, Feb 2020


Nmax = Pmax;
A = Create_Dictionary(Nmax,size(tfr,1),'Ramanujan');

tpr = zeros(size(tfr,1),size(tfr,2));

for l = 2:size(tfr,2)
    x = tfr(:,l);
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
    tpr(2:Nmax+1,l) = s;
    
end

end

