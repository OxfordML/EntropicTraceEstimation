function [ loglik ] = loglikMatrixInv_( z, invcovmat )

% computes the Gaussian loglikelihood from an observation
% vector 'z' and an inverse covariance matrix 'invcovmat'

% I would like to add a profile likelihood version of this

% take Cholesky and make sure it is positive definite
N = length(z);

logdet = -maxent_logdet(invcovmat);
quadform = z'*invcovmat*z;

loglik = -1/2*( logdet + quadform + N*log(2*pi));

end
