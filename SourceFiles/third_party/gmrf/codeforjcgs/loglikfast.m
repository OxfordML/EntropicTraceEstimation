function [ loglik,beta,tausq,logdet,quadform] = ...
    loglikfast(zarray,observedmap,i1,i2array,D,beta,parms,nu,J,profbeta,proftausq)

% computes the exact loglikelihood described in the paper
% uses symmetric matrix methods
% uses different methods for 'no nugget' and 'nugget cases'


% Name the parameters
tausq = parms(1);
kappasq = parms(2);
nug = parms(3);            % \sigma^2 = nug/tausq;
nugtrue = nug>0;


% get the grid size
boldn = size(zarray);
n1 = boldn(1);
n2 = boldn(2);
n = n1*n2;


% associate a unique index to each grid point
indsarray = reshape( 1:n, n1, n2 );


% the indices where we have observations
observedinds = indsarray( observedmap==1 );
N = length(observedinds); % number of observed values


% the number of neighbors is determined by nu
nneighbors = sum(2*(nu+1)+1-2*abs(-nu-1:nu+1));


if ~nugtrue, % no nugget case

    [Q0,neighborcount] = getQ_none(observedmap,i1,i2array,parms,nu);
    inds1 = observedinds( neighborcount(:) < nneighbors );
    inds2 = observedinds( neighborcount(:) == nneighbors );
    
    Q22 = Q0(inds2,inds2);
    cQ22 = chol(Q22,'lower');
    
    [K11,fJ] = getK11(boldn,inds1,parms,nu,J);
    

    [cK11,p] = chol(K11,'lower'); % really should add this below as well
    if p > 0,                     % for the nugget case
        loglik = -999999999;
        return;
    end
        
    
    % determinant
    logdet = 2*sum(log(diag(cK11))) - 2*full(sum(log(diag(cQ22))));
    
    Y = zeros(n,1);
    Y(observedinds) = zarray(observedinds);
    
    
    % profile out beta if requested
    if profbeta,
        M = halfsolve_nonug(D,boldn,cK11,cQ22,fJ,inds1,inds2);
        R = fullsolve_nonug(Y,boldn,cK11,Q22,fJ,inds1,inds2);
        beta = ( M'*M ) \ D'*R;
        mu = D*beta;
    else
        mu = D*beta;
    end
    
    

    
    % quadratic form
    F = halfsolve_nonug(Y-mu,boldn,cK11,cQ22,fJ,inds1,inds2);
    quadform = F'*F;
    
    
    % profile out tausq if requested    
    if proftausq,
        tausq = N/quadform;
        loglik = -1/2*logdet + N/2*log(tausq) - N/2 - N/2*log(2*pi);
    else
        loglik = -1/2*logdet -1/2*quadform - N/2*log(2*pi);
    end
    

         
else  % use a different set of calculations if there is a nugget
    
    
    % calculate the precision matrix
    Q = getQ_exact(observedmap,i1,i2array,[tausq kappasq],nu,J);


    size(Q(observedinds,observedinds))
    nnz(Q(observedinds,observedinds))
    tic;
    [cQ,~,~] = chol(Q(observedinds,observedinds),'lower','vector');
    round(toc*1000)
    A = nug/tausq*Q + speye(n);
    size(A)
    nnz(A)
    tic
    [cA,~,perm] = chol(A,'lower','vector');
    round(toc*1000)
    logdet = -2*full(sum(log(diag(cQ)))) + 2*full(sum(log(diag(cA))));

    Y = zeros(n,1);
    Y(observedinds) = zarray(observedinds);

    
    if profbeta, 
        M =   cA' \ ( cA \ ( Q*D ) );
        R = cA' \ ( cA \ ( Q*Y ) );
        beta = ( D'*M ) \ D'*R;
        mu = D*beta;
    else
        mu = D*beta;
    end

    
    W = Q*(Y-mu);
    V = cA' \ (cA \ (Y(perm)-mu(perm)) ); 
    quadform = W(perm)'*V;

    if proftausq,
        tausq = N/quadform;
        loglik = -1/2*logdet + N/2*log(tausq) - N/2 - N/2*log(2*pi);
    else
        loglik = -1/2*logdet -1/2*quadform - N/2*log(2*pi);
    end
    
end


end

