function [ loglik ] = loglikQ( zarray, observedmap,i1,i2array,D,beta, parms, nu, getQ )

% computes an approximate likelihood with approximate
% precision matrix Q

mu = D*beta;
% Makes use of loglikmatrixinv()
Q = getQ(observedmap,i1,i2array, parms, nu);
inds = observedmap==1;
loglik = loglikMatrixInv(zarray(inds)-mu(inds),Q(inds,inds));

end

