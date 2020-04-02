function [energy_s] = PD_Lasso(x,Pmax,method,lambda,n)
% Code for Ramanujan-based l1 penalized linear regression periodicity transform for 1-dim signal
% INPUT:
%  x - input signal
%  Pmax - maximal estimated period
%  method = 'Ramanujan'
%  lambda - penalty in the program
%  n - penalty function \zeta(p) = p^n
% OUTPUT:
%  energy_s - energy of period

% Written by Ziyu Chen, 2019




if(size(x,2)>size(x,1))
    x = x';
end

N = length(x);

A = Create_Dictionary(Pmax,N,method);


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


end

