function [ loglik,logdet,quadform ] = loglikK( zarray, observedmap, parms, nu, J )


tausq = parms(1);
nug = parms(3);

boldn = size(zarray);

Karray = getKarray(boldn,parms,nu,J);

Kfull = bltoeplitz2(Karray);

inds = observedmap(:)==1;
N = sum(observedmap(:));

K = Kfull(inds==1,inds==1) + nug/tausq*eye(length(inds));
cK = chol(K,'lower');

logdet = 2*sum(log(diag(cK)));
z0 = cK \ zarray(:);
quadform = z0'*z0;

loglik = -1/2*( logdet + quadform + N*log(2*pi) );

end

