function [tpr] = RPT_time_period(y,Pmax,method,lambda,nei,hop,n)
% Code for short-time Ramanujan-based l1 penalized linear regression periodicity transform for 1-dim signal
% INPUT:
%  x - input signal
%  Pmax - maximal estimated period
%  method = 'Ramanujan'
%  lambda - penalty in the program
%  nei - number of neighbors on one side
%  hop - hop in time
%  n - penalty function \zeta(p) = p^n
% OUTPUT:
%  tpr - short-time time-period representation

% Written by Ziyu Chen, March 2020




if(size(y,2)>size(y,1))
    y = y';
end
N = length(y);
seq = 1:hop:N;
le = length(seq);
tpr = zeros(Pmax,le);

Nmax = Pmax;
fprintf(['Columns total: ',num2str(le),'; now:     ']) ;
for l = 1:le
    fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',l) ; fprintf([tmp]) ;
    idx = max(1,seq(l)-nei):min(seq(l)+nei,N);
    L = length(idx);
    x = y(idx);
    A = Create_Dictionary(Nmax,L,method);
    
    
    % for i=1:size(A,2)
    %     A(:,i) = A(:,i)./norm(A(:,i));
    % end
    
    % Penalty Vector Calculation
    
    penalty_vector = [];
    for i=1:Nmax
        k = 1:i;
        k_p = k(gcd(k,i)==1);
        k_p = length(k_p);
        penalty_vector=cat(1,penalty_vector,i*ones(k_p,1));
    end
    
    penalty_vector = penalty_vector.^(n);
    
    D = zeros(length(penalty_vector));
    DD = zeros(length(penalty_vector));
    
    for i = 1:length(penalty_vector)
        D(i,i) = penalty_vector(i);
        DD(i,i) = 1/penalty_vector(i);
    end
    
    s = lasso(sqrt(N)*A*DD,sqrt(N)*x,'Lambda',lambda,'Standardize',false);
    %s = D*s;
    
    energy_s = 0.*[1:Nmax];
    current_index_end = 0;
    for i=1:Nmax
        i_all = 1:i;
        m = i_all(gcd(i_all,i)==1);
        current_index_start = current_index_end + 1;
        current_index_end = current_index_end + size(m,2);
        
        for j=current_index_start:current_index_end
            
            energy_s(i) = energy_s(i)+((abs(s(j)))^2);
            %energy_s(i) = energy_s(i)+abs(s(j));
        end
        
    end
    
    tpr(:,l) = energy_s;
    
end
end

