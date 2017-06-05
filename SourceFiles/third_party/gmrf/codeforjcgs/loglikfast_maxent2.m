function loglik = ...
    loglikfast_maxent2(zarray,observedmap,i1,i2array,D,beta,parms,nu,J,profbeta,proftausq)

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


% calculate the precision matrix
Q = getQ_exact(observedmap,i1,i2array,[tausq kappasq],nu,J);

B = Q(observedinds,observedinds);
tic;
l1 = maxent_logdet(B)
%[cQ,~,~] = chol(B,'lower','vector');
round(toc*1000)

A = Q(observedinds,observedinds) + speye(length(observedinds));
tic;
l2 = log(nug/tausq)*size(A) + maxent_logdet(A)
%[cA,~,perm] = chol(A,'lower','vector');
round(toc*1000)

%log1 = l2
%log1t = 2*full(sum(log(diag(cA))))
%diff1 = abs(l1 - 2*full(sum(log(diag(cQ))))) / abs(2*full(sum(log(diag(cQ)))))
%diff2 = abs(l2 - 2*full(sum(log(diag(cA))))) / abs(2*full(sum(log(diag(cA)))))
logdet = l2 + l1;%-2*full(sum(log(diag(cQ)))) + 2*full(sum(log(diag(cA))));

Y = zeros(n,1);
Y(observedinds) = zarray(observedinds);


% if profbeta, 
%     M =   cA' \ ( cA \ ( Q*D ) );
%     R = cA' \ ( cA \ ( Q*Y ) );
%     beta = ( D'*M ) \ D'*R;
%     mu = D*beta;
% else
%     mu = D*beta;
% end


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

