function [i1,i2array] = getSparsePattern_periodic(observedmap,nu)

% get size of the grid and missingness pattern
[n1,n2] = size(observedmap);
n = n1*n2;
indsarray = reshape( 1:n, n1, n2 );
observedinds = indsarray( observedmap==1 );
N = length(observedinds);

% "width" of the neighborhood
len = 2*(nu+1)+1;

% the values of the paramters are not important here
% we just use them to generate the neighborhood pattern
tau = 1;
lam = 0.99;
basecoefs = zeros(len);
basecoefs(1,1) = tau;
basecoefs([2 len len+1 len*(len-1)+1]) = -tau*lam/4;
coefmat = basecoefs;
% convolution
if nu >=1
    for j = 1:nu
        coefmat = real(ifft2( fft2(basecoefs).*fft2(coefmat) ));
    end
end
% These are the neighbors with non-zero coefficients
neighborhood = abs(coefmat) > 1e-12; 
nneighbors = sum(neighborhood(:));

% each observation has potentially nneighbors-1 neighbors,
% and the precision matrix also has a non-zero entry
% on the diagonal, so each row of Q has nneighbors nonzero entries

% i1 stores the nonzero row indices of Q
i1 = reshape( repmat( observedinds', nneighbors,1 ) , nneighbors*N,1);
% i2 stores the nonzero column indices of Q
i2array = NaN(nneighbors,N);

% In this periodic version, we allow grid points on the opposite
% side of the grid to be neighbors. We still have to use the 'sh'
% mechanism to shift neighborinds and neighbormap to 
% match neighborhood. But the code is simpler here than
% it is in getSparsePattern()
for j = 1:length(observedinds),
    k = observedinds(j);
    j2 = ceil( k/n1 );
    j1 = k - (j2-1)*n1;
    r1 = mod((j1-nu-1:j1+nu+1)-1,n1)+1;
    r2 = mod((j2-nu-1:j2+nu+1)-1,n2)+1;
    neighborinds =  indsarray(r1,r2);
    neighbormap = observedmap(r1,r2);
    sh = [nu+2:2*(nu+1)+1,1:nu+1];
    shiftinds = neighborinds(sh,sh);
    shiftmap = neighbormap(sh,sh);
    i2array(:,j) = shiftinds(neighborhood==1).*shiftmap(neighborhood==1);
end



end