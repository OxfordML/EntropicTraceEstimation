function [ zarray ] = simulateMRF( boldn, parms, nu, J )

% Simulates a Markov Random Field on a grid of size boldn

tausq = parms(1);
kappasq = parms(2);
nug = parms(3);

n1 = boldn(1);
n2 = boldn(2);

% define the frequency arrays
w1vec = (0:n1*J-1)'*2*pi/(n1*J);
w2vec = (0:n2*J-1)'*2*pi/(n2*J);
w1 = repmat( w1vec, 1, n2*J );
w2 = repmat( w2vec, 1, n1*J )';

% compute the spectral density
fJ = 1/tausq*( nug + (kappasq+4-2*( cos(w1) + cos(w2) ) ).^(-nu-1) );

% use FFTs to simulate periodic process on a grid of size
% (n1*J,n2*J)
zfull = real( ifft2( sqrt(fJ).*fft2(randn(n1*J,n2*J)) ) );
% extract the simulated values on the (n1,n2) grid
zarray = zfull(1:n1,1:n2);

end

