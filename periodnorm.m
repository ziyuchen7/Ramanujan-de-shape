function y=periodnorm(x)

% Calculates the norm of a vector using the inner product
% defined on the p-periodic subspaces:
%
% ||x|| = sqrt(<x,x>) = norm(x)/sqrt(length(x))
%
% syntax: y = periodnorm(x)
% input:  x = vector to find norm of
% output: y = (p-periodic) norm of x 
%
% See Sethares and Staley, "The Periodicity Transform"
% IEEE Trans. Signal Processing, 1998.

y=norm(x)/sqrt(length(x));
