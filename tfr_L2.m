function [tpr] = tfr_L2(tfr,Pmax,method)
% Code for TFR using constrained L1 minimization
% INPUT:
%  tfr - input time frequency representation
%  Pmax - bin index of predetermined largest frequency
%  method = 'Ramanujan'
% OUTPUT:
%  tpr - result TPR
% Here the penalty matrix D is chosen as f(P_i) = (P_i)^2 - line 41

% Written by Ziyu Chen, 2019 based on
% Strength_vs_Period_L2(x,Pmax,method) by S.V. Tenneti and P. P. Vaidyanathan

Nmax = Pmax;
A = Create_Dictionary(Nmax,size(tfr,1),method);

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
    for i=1:Nmax
        k=1:i;k_red=k(gcd(k,i)==1);k_red=length(k_red);
        penalty_vector=cat(1,penalty_vector,i*ones(k_red,1));
    end
    
    penalty_vector = penalty_vector.^2;
    
    D = diag((1./penalty_vector).^2);
%     B = A*D*A';
%     [U,S,V] = svd(B);
%     invS = zeros(size(x,1),size(x,1));
%     for i=1:size(x,1)
%         if S(i,i)>1e-6
%             invS(i,i) = 1/S(i,i);
%         else
%             invS(i,i) = 0;
%         end
%     end
%     
%     invB = V*invS*U';
%     PP = D*A'*invB;
    PP = D*A'*inv(A*D*A');
    s = PP*x;
    %s = DD*s;
    
    energy_s = 0.*[1:Nmax];
    current_index_end = 0;
    for i=1:Nmax
        k_orig = 1:i;k=k_orig(gcd(k_orig,i)==1);
        current_index_start = current_index_end + 1;
        current_index_end = current_index_end + size(k,2);
        
        for j=current_index_start:current_index_end
            
            energy_s(i) = energy_s(i)+((abs(s(j)))^2);
            %energy_s(i) = energy_s(i)+abs(s(j));
        end
        
    end
    
    %energy_s(1) = 0;
    
    tpr(2:Nmax+1,l) = energy_s;
end

end
