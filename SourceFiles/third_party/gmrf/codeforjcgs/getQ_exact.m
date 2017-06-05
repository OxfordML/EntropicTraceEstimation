function [ Q,K11,cK11,cQ22,neighborcount ] = getQ_exact(observedmap,i1,i2array,parms,nu,J)

% returns the exact precision matrix for Z

tausq = parms(1);
kappasq = parms(2);

% get the grid size
boldn = size(observedmap);
n1 = boldn(1);
n2 = boldn(2);
n = n1*n2;

% associate a unique index to each grid point
indsarray = reshape( 1:n, n1, n2 );

% the indices where we have observations
observedinds = indsarray( observedmap==1 );

% the number of neighbors is determined by nu
nneighbors = sum(2*(nu+1)+1-2*abs(-nu-1:nu+1));

[Q0,neighborcount] = getQ_none(observedmap,i1,i2array,[tausq kappasq],nu);

% indicate the partially and fully neighbored observations
inds1 = observedinds( neighborcount(:) < nneighbors );
inds2 = observedinds( neighborcount(:) == nneighbors );

Q22 = Q0(inds2,inds2);
Q21 = Q0(inds2,inds1);


K11 = getK11(boldn,inds1,parms,nu,J);

cK11 = chol(K11,'lower');
Kinv = cK11' \ ( cK11 \ speye(length(inds1)) );

%B = cK11 \ speye(length(inds1));
%Kinv = B'*B;

[cQ22,~,perm] = chol(Q22,'lower','vector');


R1 = cQ22 \ Q21(perm,:);
R2 = R1'*R1;
Q11 = Kinv + R2;

Q = Q0;
Q(inds1,inds1) = Q11;



end

