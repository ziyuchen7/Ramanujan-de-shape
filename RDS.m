function [tpr] = RDS(tfr,Pmax,method,lambda,n)
% Code for Ramanujan de-shape (RDS)
% INPUT:
%  tfr - input time frequency representation
%  Pmax - bin index of predetermined largest frequency
%  method = 'Ramanujan'
%  lambda - penalty in the program 
%  n - penalty function \zeta(p) = p^n
% OUTPUT:
%  tpr - RDS

% Written by Ziyu Chen, December 2019

A = Create_Dictionary(Pmax,size(tfr,1),method);

tpr = zeros(size(tfr,1),size(tfr,2));
N = size(tfr,1);

fprintf(['Columns total: ',num2str(size(tfr,2)),'; now:     ']) ;

for l = 1:size(tfr,2)
    fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',l) ; fprintf([tmp]) ;
    x = tfr(:,l);
    
    
    
    
    % for i=1:size(A,2)
    %     A(:,i) = A(:,i)./norm(A(:,i));
    % end
    
    % Penalty Vector Calculation
    
    penalty_vector = [];
    for i = 1:Pmax
        k = 1:i;
        k_p = k(gcd(k,i)==1);
        k_p = length(k_p);
        penalty_vector = cat(1,penalty_vector,i*ones(k_p,1));
    end
    
    penalty_vector = penalty_vector.^(n);
    
    D = zeros(length(penalty_vector));
    DD = zeros(length(penalty_vector));
    
    for i = 1:length(penalty_vector)
        D(i,i) = penalty_vector(i);
        DD(i,i) = 1/penalty_vector(i);
    end
    
    s = lasso(sqrt(N)*A*DD,sqrt(N)*x,'Lambda',lambda,'Standardize',false);
    
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

