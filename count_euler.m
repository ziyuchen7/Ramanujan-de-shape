function [s] = count_euler(N)
% Totient summation function \Phi(N)
% By Ziyu Chen, 2019

if N==0
    s = 0;
    return 
end

s = 0;
for i = 1:N
    k = 1:i;
    sub_k = k(gcd(k,i)==1);
    l = length(sub_k);
    s = s+l;
end

end


