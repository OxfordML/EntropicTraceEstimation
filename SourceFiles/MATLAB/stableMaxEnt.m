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

for k = 1:100000 %arbitrary loop size
    p_tilde = exp(A*n-1);
    
    i = mod(k, size(mu,1)) + 1;
    
    lambda = log(mu(i)/(A(:,i)'*p_tilde));
    
    n(i) = n(i) + lambda;
    
end

a = repmat(exp(A*n - 1), 1, size(x,1));
moments = abs(sum(A'.*a,1))';
p = exp(A*n-1);
end

%shallow_water1 - 10 moments
%Truth = -5.7151e+04;
%Est = -5.7322e+04;
%rel error: 0.003

%shallow_water2 - 10 moments
%Truth =  -9.1731e+04
%Est = -9.2203e+04
%rel error: 0.0051

%apache1 - 10 moments
%Truth = -3.2090e+05;
%Est = -3.1908e+05;
%rel error: 0.0057

%finan512 - 10 moments
%Truth = -2.0942e+05;
%Est = -2.1300e+05;
%rel error: 0.0171

%obstclae - 10 moments
%Truth = -3.5186e+04
%Est = -3.5851e+04
%rel error: 0.0026 

%jnlbrng1
%Truth = -5.3877e+04
%Est = -5.5001e+04
%rel error: 0.0158

%$ 100^2 $ &   0.05 &   0.14 &   0.00 &   0.17 &   0.21 &   0.00 \\ 
%$ 200^2 $ &   0.34 &   0.38 &   0.00 &   0.68 &   0.67 &   0.00 \\ 
%$ 300^2 $ &   1.27 &   0.76 &   0.00 &   2.33 &   1.44 &   0.00 \\ 
%$ 400^2 $ &   4.62 &   1.52 &   0.00 &   9.08 &   2.51 &   0.00 \\ 
%$ 500^2 $ &   9.46 &   2.06 &   0.00 &  28.43 &   4.05 &   0.00 \\ 
%$ 600^2 $ &  24.51 &   3.26 &   0.00 & 306.23 &   6.69 &   0.00 \\ 
%$ 700^2 $ & 135.71 &   4.53 &   0.00 & 517.71 &   8.32 &   0.00 \\ 
%$ 800^2 $ & 349.09 &   9.91 &   0.00 & 593.88 &  12.21 &   0.00 \\ 
%$ 900^2 $ & 478.09 &  12.54 &   0.00 & 3282.25 &  17.07 &   0.00 \\ 
%$ 1000^2 $ & 973.52 &  11.08 &   0.00 & 2212.15 &  17.98 &   0.00 \\

