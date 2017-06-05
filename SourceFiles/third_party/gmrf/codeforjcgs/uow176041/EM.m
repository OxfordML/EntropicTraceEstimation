% EM algorithm to find MLEs of K and sigma2 (fine-scale variance)

function [K_em sig2_em T diff]=EM(S,z,V,sigeps,V2,maxit,avgtol)

% V2 is the diagonal covariance matrix for the fine-scale variation
%     (default: V2=eye(size(S,1));  )

n=size(S,1);

% default values for the last two parameters
if nargin<7, avgtol=1e-6; end
if nargin<6, maxit=200; end
if nargin<5, V2=sparse(1:n,1:n,1); end

diagV=diag(V);
diagV2=diag(V2);

% initial values
varest=var(z,1);
K_old=.9*varest*eye(size(S,2));
sig2=.1*varest;
t=1;
done=0;
    
while done==0,
   
    % update help terms
    diagDinv=(sig2(t)*diagV2+sigeps*diagV).^(-1);
    DInv=sparse(1:n,1:n,diagDinv);
    tempt=inv(inv(K_old)+S'*DInv*S);
    
    % update K
    SigInv2=(tempt*S')*DInv;
    KSDInv=K_old*S'*DInv;
    KSSigInv=KSDInv-KSDInv*S*SigInv2;
    muEta=KSSigInv*z;
    SigEta=K_old-KSSigInv*S*K_old;
    K_new=SigEta+muEta*muEta';

    % update sigma_xi (sig2)
    muEps=sig2(t)*(DInv*z-DInv*S*(SigInv2*z));
    trSigInv=trace(DInv)-trace(SigInv2*DInv*S);
    sig2(t+1)=1/n*(n*sig2(t)-(sig2(t))^2*trSigInv+muEps'*muEps);
    
    % check for convergence
    diff=sum(sum((K_new-K_old).^2,1),2)+(sig2(t+1)-sig2(t))^2;
    if diff<avgtol*(size(S,2))^2, done=1; end
    if t>maxit, 
        done=1; 
        disp(strcat('Algorithm did not converge after ',num2str(maxit),' iterations')); 
    end
    
    t=t+1;
    K_old=K_new;

end

K_em=K_new;
sig2_em=sig2(t);
T=t;