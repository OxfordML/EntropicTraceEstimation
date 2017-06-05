function [ K11, fJ ] = getK11( boldn, inds1, parms, nu, J )

% We have to do this in multiple parts of the code

n1 = boldn(1);
n2 = boldn(2);
n = n1*n2;

% associate a unique index to each grid point
indsarray = reshape( 1:n, n1, n2 );


tausq = parms(1);
kappasq = parms(2);

% compute the array of covariances
[Karray,fJ] = getKarray(boldn,[tausq kappasq],nu,J);

% Kmat will be the appropriate columns
% of the covariance matrix
K11 = zeros(length(inds1),length(inds1));
% loop over the columns of K11
for obs = 1:length(inds1)
    % create vector containing covariances of observation 
    % inds1(obs) with all other obs.
    col = inds1(obs);
    % get the coordinates of the grid locations (p1,p2)
    p2 = ceil( col/n1 );
    p1 = col - (p2-1)*n1;
    % identify the indices of Karray to extract for
    % constructing the column of the covariance matrix
    q1 = abs( (-p1+1 : -p1+n1)' ) + 1;
    q2 = abs( (-p2+1 : -p2+n2)' ) + 1;
    ilist = indsarray(q1,q2);
    % assign the 'obs' column of Kmat
    K11(:,obs) = Karray(ilist(inds1));
end

end

