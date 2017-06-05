% Timing experiment to compare approximate and
% exact likelihoods



% grid sizes and nu parameters
nvec = 500;
nuvec = 0;

% amount of oversampling to compute covariances
J = 2;

% array to save timing results
timemat = zeros(length(nvec),3,length(nuvec));

% loop over values of nu
nu = nuvec;

% we compute three likelihoods for each (parameter,gridsize)
% combination--one with the approximation and no adjustment
% to the precision matrix, one with the approximation 
% and a periodic adjustment to the precision matrix
% and one exact likelihood.
getQ = @getQ_none;
pattern1 = @getSparsePattern;
pattern2 = @getSparsePattern_periodic;
pattern3 = @getSparsePattern;


% loop over grid sizes
n1 = nvec;
n2 = n1;
n = n1*n2;
boldn = [n1,n2];
observedmap = zeros(n1,n2);
observedmap(:) = 1;
zarray = simulateMRF(boldn,parms,nu,2*J);
D = ones(n,1);
beta = 0;

[i11,i21array] = pattern1(observedmap,nu);
[i12,i22array] = pattern2(observedmap,nu);
[i13,i23array] = pattern3(observedmap,nu);
i1 = i13;
i2array = i23array;

for j=1:2
    for i = 1:12
        % parameters to use
        tausq = 1;
        kappasq = 1/10;
        nug = 0.0;
       
        parms = [tausq kappasq nug];
        parms(j) = i/12;

        lq(i, j) = loglikQ(zarray,observedmap,i11,i21array,D,beta,parms(1:2),nu,getQ);
        lq_(i, j) = loglikQ_(zarray,observedmap,i11,i21array,D,beta,parms(1:2),nu,getQ);
    end
end


% print out the results

