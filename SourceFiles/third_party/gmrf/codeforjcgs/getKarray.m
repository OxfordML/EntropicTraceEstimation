function [ Karray, fJ ] = getKarray( boldn, parms, nu, J )

% Returns an array of covariances K[i,j] = K((i,j))
% on a grid of size boldn

% We use an approximation that relies on computing a 
% 2d DFT on grid that is a factor of J larger in both dimensions

% parameters
tausq = parms(1);
kappasq = parms(2);

% grid dimensions
n1 = boldn(1);
n2 = boldn(2);

% arrays of frequencies
w1vec = (0:n1*J-1)'*2*pi/(n1*J);
w2vec = (0:n2*J-1)'*2*pi/(n2*J);
cosw1vec = cos(w1vec);
cosw2vec = cos(w2vec);
cosw1 = repmat( cosw1vec, 1, n2*J );
cosw2 = repmat( cosw2vec, 1, n1*J )';

% spectral density
fJ = 1./(kappasq + 4 - 2*( cosw1 + cosw2 ) );
if nu>=1 
    for j = 1:nu,
        fJ = fJ.*fJ;
    end
end
fJ = 1/tausq*(fJ);


% take inverse DFT to get covariances
KarrayFull = real(ifft2(fJ));

% extract covariances on the grid we need
Karray = KarrayFull(1:n1,1:n2);

end

