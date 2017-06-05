function [p, n, moments] = stableMaxEnt(mu, x)

%%
% stableMaxEnt(mu, x)
% mu is the moment information
% x are the sample points used 
%

n = zeros(size(mu));

A = ones(length(x), size(mu,1));
for i = 2 : size(mu,1)
   A(:, i) = (x.^(i - 1));
end

for k = 1:10000 %arbitrary loop size
    p_tilde = exp(A*n-1);
    
    i = mod(k, size(mu,1)) + 1;
    
    lambda = log(mu(i)/(A(:,i)'*p_tilde));
    
    n(i) = n(i) + lambda;
    
end

%a = repmat(exp(A*n - 1), 1, size(x,1));
%moments = abs(sum(A'.*a,1))';
p = exp(A*n-1);
moments = 0;
end