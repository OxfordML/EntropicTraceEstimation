function [ loglik,tausq ] = loglikvecchia( z, coords, blocks,D0,beta, parms, nu, J )
% computes the block composite likelihood for MRF model with nugget

N = length(z);
% coords must have integer values

tausq = parms(1);
kappasq = parms(2);
nug = parms(3);



mini = min(coords(:,1));
maxi = max(coords(:,1));
minj = min(coords(:,2));
maxj = max(coords(:,2));

% form the covariance array
boldn = [maxi-mini+1,maxj-minj+1];
Karray = getKarray(boldn,[tausq kappasq],nu,J);
Karray(1,1) = Karray(1,1)+nug/tausq;

% compute the likelihood on each block
nblocks = length(blocks);
blogdet = zeros(nblocks,1);
bquadform = zeros(nblocks,1);
z0 = z - D0*beta;


% exact likelihood for first block
b = 1;
binds = blocks{b};
minib = min(coords(binds,1));
maxib = max(coords(binds,1));
minjb = min(coords(binds,2));
maxjb = max(coords(binds,2));
bcoords = coords(blocks{b},:);
bmap = zeros(maxib-minib+1,maxjb-minjb+1);
for j = 1:length(binds)
    bmap( bcoords(j,1)-minib+1,bcoords(j,2)-minjb+1 ) = j;
end
bKarray = Karray( 1:(maxib-minib+1), 1:(maxjb-minjb+1) );
bcovmatfull = bltoeplitz2( bKarray );
bcovmat = bcovmatfull( bmap>0, bmap>0 );
sinds = bmap( bmap > 0 );
[~,sorder] = sort(sinds);
cbcovmat = chol( bcovmat(sorder,sorder), 'lower' );
z2 = z0(binds);
bzstar = cbcovmat \ z2;
blogdet(b) = 2*sum(log(diag(cbcovmat)));
bquadform(b) = (bzstar'*bzstar);


for b = 2:nblocks
    
    N1 = length(blocks{b-1});
    N2 = length(blocks{b});
    binds = [blocks{b-1};blocks{b}];
    minib = min(coords(binds,1));
    maxib = max(coords(binds,1));
    minjb = min(coords(binds,2));
    maxjb = max(coords(binds,2));
    bcoords = coords(binds,:);
    
    bmap = zeros(maxib-minib+1,maxjb-minjb+1);
    for j = 1:length(binds)
        bmap( bcoords(j,1)-minib+1,bcoords(j,2)-minjb+1 ) = j;
    end

    bKarray = Karray( 1:(maxib-minib+1), 1:(maxjb-minjb+1) );
    bcovmatfull = bltoeplitz2( bKarray );
    bcovmat = bcovmatfull( bmap>0, bmap>0 );
    sinds = bmap( bmap > 0 );
    [~,sorder] = sort(sinds);
    cbcovmat = chol( bcovmat(sorder,sorder), 'lower' );
    z2 = z0(binds);

    bzstar = cbcovmat \ z2;
    blogdet(b) = 2*sum(log(diag(cbcovmat(N1+1:N1+N2,N1+1:N1+N2))));
    bquadform(b) = (bzstar(N1+1:N1+N2)'*bzstar(N1+1:N1+N2));
    
end

tausq = N/sum(bquadform);
loglik = -1/2*sum(blogdet) + N/2*log(tausq) - N/2 - N/2*log(2*pi);


end