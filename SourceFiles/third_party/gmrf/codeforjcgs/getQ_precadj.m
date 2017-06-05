function [ Q,neighborcount ] = getQ_precadj(observedmap,i1,i2array,parms,nu)

% returns the approximate precision matrix with an adjustment to
% the precision of non-fully neighbored observations

% this should only be used if the inverse covariances are
% diagonally dominant, as they are in the nu=0 case

% define the parameters
tausq = parms(1);
kappasq = parms(2);

% get the grid sizes and observed pattern
[n1,n2] = size(observedmap);
n = n1*n2;
indsarray = reshape( 1:n, n1, n2 );
observedinds = indsarray( observedmap==1 );
N = length(observedinds);

% This is the "width" of the neighborhood
len = 2*(nu+1)+1;

% the "theta" for the nu=0 case are given by basecoefs
% the (1,1) entry of basecoefs is for h=(0,0)
% the neighborhood here is periodic, so the (len,1) entry
% corresponds to h=(-1,0)
basecoefs = zeros(len);
basecoefs(1,1) = (kappasq+4);
basecoefs([2 len len+1 len*(len-1)+1]) = -1;

% if nu >= 1 , the coefficients are convolutions
% of basecoefs with itself.
% this is why we use the periodic ordering in basecoefs
coefmat = basecoefs;
if nu >=1
    for j = 1:nu
        coefmat = real(ifft2( fft2(basecoefs).*fft2(coefmat) ));
    end
end
% these are the values of h for which the coefficient is nonzero
neighborhood = abs(coefmat) > 1e-12;
nneighbors = sum(neighborhood(:));
% finally, multiply by tausq
coefvec = tausq*coefmat(neighborhood==1);

% this stores all of the entries of q. The jth column of 
% Qarray stores all of the entries that go in the jth column of Q
Qarray = repmat( coefvec, 1, N );

% an=0 tells us when a neighbor is not in the list of locations
an = i2array ~= 0;
% count up the number of available neighbors for each observation
neighborcount = sum(an)';
% extract only the entries for which we have a neighbor
% and vectorize
i1 = i1(an);
i2 = i2array(an);
% precision adjustment
Qarray(1,:) = Qarray(1,:).*( (neighborcount-1)'/(nneighbors-1) );
Qvals = Qarray(an);

% assign to Q
Q = sparse(i1,i2,Qvals,n,n);


end

