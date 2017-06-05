% This function does the binned frobenius norm estimation for K

function K=binest(data,S,sige,V,sigxi,V2)

n=size(S,1);
if nargin<6, V2=sparse(1:n,1:n,1); end

global SigHat D_ eigs Us lam0;

lon=data(:,2); lat=data(:,1); res=data(:,3); 


% bin the data locations

W=bin(lon,lat); 
% (replace this with your own binning procedure for use on a 
%   different dataset)


% calculate the binned quantities

M=length(W);

for j=1:M,
    
    wj=zeros(n,1); wj(W{j})=1;
    nj=sum(wj);
    
    S_(j,:)=wj'*S/nj;
    V_(j,j)=wj'*V*wj/nj;
    V2_(j,j)=wj'*V2*wj/nj;
    
    Z_(j,1)= mean(res(W{j}));
    Z2_(j,1)= mean((res(W{j})).^2);
    A_(j,j)= sqrt(nj)/Z2_(j,1);
    
end

for i=1:M,
    for j=1:M,
        
        if i==j, SigHat(i,j)=Z2_(i,1) ; 
        else SigHat(i,j)= Z_(i,1)*Z_(j,1);
        end
        
    end
end

D_=sigxi*V2_+sige*V_;


%%% lifting

A=D_^(-.5)*(SigHat-D_)*D_^(-.5);
[U,Lam] = eig(A);
[eigs ind]=sort(diag(Lam),'descend');
Us=U(:,ind);

if min(eigs)<0,
    
    q=size(S,2); % q = r;
    a=-1;
    
    while a<0,
        
        disp(strcat('trying q=',num2str(q)))
        if q>=0,
            lam0=max([quantile(eigs,(M-q)/M) 0]); 
        else lam0=eigs(1)-q;
        end
        
        try
            a=fzero(@trvar,[0 10]);
        catch exception
            q=q-1;            
        end  
            
    end

    eigs_star=eigs;
    eigs_star(eigs<lam0)=lam0.*exp(a.*(eigs(eigs<lam0)-lam0));
    SigHat_star= D_^(.5)* Us*diag(eigs_star)*Us' *D_^(.5)+D_ ;

else 
    SigHat_star=SigHat;
    
end


% QR decomposition and estimation of K

[Q,R] = qr(S_,0);
K=(R\Q')*(SigHat_star-D_)*(Q/R');
