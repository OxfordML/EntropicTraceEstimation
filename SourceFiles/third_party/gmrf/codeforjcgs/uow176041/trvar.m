%%  function to find a which preserves the total variability

function diffvar=trvar(a)

global SigHat D_ eigs Us lam0;

eigs_star=eigs;
eigs_star(eigs<lam0)=lam0.*exp(a.*(eigs(eigs<lam0)-lam0));

A_star=Us*diag(eigs_star)*Us';

SigHat_star=D_^(.5)*A_star*D_^(.5)+D_ ;

diffvar=(trace(SigHat_star-SigHat));

