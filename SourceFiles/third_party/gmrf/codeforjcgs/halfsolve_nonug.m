function [ outmat ] = halfsolve_nonug( X, boldn, cK11, cQ22, fJ, inds1, inds2  )

% if LL^T = \Sigma (no nugget), we compute here
%
% L^{-1}X 
%
% for an nxp matrix X

n1 = boldn(1);
n2 = boldn(2);
n = n1*n2;


[n0,p] = size(X);
if n0 ~= n, error('length of X not equal to number of grid points'); end

outmat = zeros(size(X));
for j = 1:p
    V = zeros(n,1);
    V(inds1) = cK11' \ (cK11 \ X(inds1,j) );
    Varray = zeros(size(fJ));
    Varray(1:n1,1:n2) = reshape(V,n1,n2);
    Warray = real(ifftn( fJ.*fftn(Varray) ));
    W0array = Warray(1:n1,1:n2);
    
    E = zeros(n,1);
    E(inds1) = X(inds1,j);
    E(inds2) = X(inds2,j) - W0array(inds2);
    
    F = zeros(n,1);
    F(inds1) = cK11 \ E(inds1);
    F(inds2) = cQ22' * E(inds2);
    
    outmat(:,j) = F;
end


end

