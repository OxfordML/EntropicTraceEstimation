function [ loglik ] = loglikMatrixInv2( z, invcovmat )

% computes the Gaussian loglikelihood from an observation
% vector 'z' and an inverse covariance matrix 'invcovmat'

% I would like to add a profile likelihood version of this

% take Cholesky and make sure it is positive definite
[cholmat,p] = chol(invcovmat,'lower');
N = length(z);
% if not PD, then return a really low loglikelihood
if p > 0,
    loglik = -9999999;
    return;
end

logdet = -2*sum(log(diag(cholmat)))
logdet_ = -maxent_logdet(invcovmat)
quadform = z'*invcovmat*z;

loglik = -1/2*( logdet + quadform + N*log(2*pi));

end
