function [tpr] = vector_RDS(tfr,Pmax,method,lambda,nei,n)
% Code for Vectorized Ramanujan de-shape
% INPUT:
%  tfr - input time frequency representation
%  Pmax - bin index of predetermined largest frequency
%  method = 'Ramanujan'
%  lambda - penalty in the program 
%  nei - number of neighbors on one side chosen at each time
%  n - penalty function \zeta(p) = p^n
% OUTPUT:
%  tpr - VRDS

% Written by Ziyu Chen, March 2020


A = Create_Dictionary(Pmax,size(tfr,1),method);

tpr = zeros(size(tfr,1),size(tfr,2));
N = size(tfr,1);

fprintf(['Columns total: ',num2str(size(tfr,2)),'; now:     ']) ;

for l = 1:size(tfr,2)
    fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',l) ; fprintf([tmp]) ;
    idx = max(1,l-nei):min(size(tfr,2),l+nei);
    y = [];
    L = length(idx);
    AA = [];
    for k = 1:L
        AA = [AA; A];
        y = [y; tfr(:,idx(k))];
    end
    
    
    % for i=1:size(A,2)
    %     A(:,i) = A(:,i)./norm(A(:,i));
    % end
    
    % Penalty Vector Calculation
    
    penalty_vector = [];
    for i=1:Pmax
        z = 1:i;
        z_p = z(gcd(z,i)==1);
        z_p = length(z_p);
        penalty_vector = cat(1,penalty_vector,i*ones(z_p,1));
    end
    
    penalty_vector = penalty_vector.^(n);
    
    D = zeros(length(penalty_vector));
    DD = zeros(length(penalty_vector));
    
    for i = 1:length(penalty_vector)
        D(i,i) = penalty_vector(i);
        DD(i,i) = 1/penalty_vector(i);
    end
    
    s = lasso(sqrt(L*N)*AA*DD,sqrt(L*N)*y,'Lambda',lambda,'Standardize',false);
    
    %s = DD*s;
    
    energy_s = 0.*[1:Pmax];
    current_index_end = 0;
    for i=1:Pmax
        i_all = 1:i;
        m = i_all(gcd(i_all,i)==1);
        current_index_start = current_index_end + 1;
        current_index_end = current_index_end + size(m,2);
        
        for j=current_index_start:current_index_end
            
            energy_s(i) = energy_s(i)+((abs(s(j)))^2);
            %energy_s(i) = energy_s(i)+abs(s(j));
        end
        
    end
    
    %energy_s(1) = 0;
    
    tpr(2:Pmax+1,l) = energy_s;
end

end

