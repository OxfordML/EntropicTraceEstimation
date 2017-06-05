% Analysis of AOT data
set(0,'defaultaxesfontsize',14)
set(0,'defaulttextfontsize',14)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'defaulttextfontname','Times New Roman')


%cd \\stat.ad.ncsu.edu\Redirect\guinness\Documents\research\carlik\code\codeforjcgs

load('redseadata.mat')


imagesc(subdata)
boldn = size(subdata);
n1 = boldn(1);
n2 = boldn(2);
n = n1*n2;
observedmap = double(observedmap);
zarray = double(subdata);


nu = 1;

J = 4;

getQ = @getQ_none;
pattern2 = @getSparsePattern;

[i1,i2array] = pattern2(observedmap,nu);

sum( sum(i2array~=0) ~= 13 )
N = sum(observedmap(:));

D = ones(n,1);
beta = 0.1;

% nugget
options = optimset('display','iter','tolFun',1e-4);
tic
%-loglikK(zarray,observedmap,[1000 exp(log([0.01 0.3]))],nu,J);
%-loglikfast_maxent(zarray,observedmap,i1,i2array,D,.1,[1000 exp(log([0.01 0.3]))],nu,J,1,1);
% funtomax = @(x) -loglikfast(zarray,observedmap,i1,i2array,D,.1,[1000 exp(x)],nu,J,1,1);
% [maxparms,~] = fminunc(funtomax,log([0.01 0.3]),options);
% kappasq = exp(maxparms(1));
% nug = exp(maxparms(2));
% [loglik0,betahat,tausq] = loglikfast(zarray,observedmap,i1,i2array,D,beta,[1 kappasq nug],nu,J,1,1);
% toc
% fprintf('%4.4f & %6.2f & %7.4f & %7.4f \n',betahat,sqrt(tausq),sqrt(kappasq),sqrt(nug/tausq))


% no nugget
% options = optimset('display','iter','tolFun',1e-4);
% tic
% funtomax = @(x) -loglikfast(zarray,observedmap,i1,i2array,D,0.1488435,[6035.303 exp(x) 0],nu,J,1,1);
% [maxparms00,~] = fminunc(funtomax,log([0.01305645]),options);
% kappasq00 = exp(maxparms00(1));
% [loglik00,betahat00,tausq00] = loglikfast(zarray,observedmap,i1,i2array,D,beta,[1 kappasq00 0],nu,J,1,1);
% toc
% fprintf('%4.4f & %6.2f & %7.4f & %7.4f \n',betahat00,sqrt(tausq00),sqrt(kappasq00),sqrt(0/tausq00))
% loglik0-loglik00
% 
% 
% % from INLA
% % max.edge=8
% betahat1 = 0.1608928;
% parms1 = [11938.83 0.001857262 1.048826];
% loglik1 = loglikfast(zarray,observedmap,i1,i2array,D,betahat1,parms1,nu,J,0,0);
% t1 = 0.61;
% llINLA1 = 69515.71;
% llINLA1 - loglik1
% fprintf('%6.1f \n',loglik0-loglik1)
% fprintf('%4.4f & %6.2f & %7.4f & %7.4f \n',betahat1,sqrt(parms1(1)),sqrt(parms1(2)),sqrt(parms1(3)/parms1(1)))
% 
% % max.edge = 6
% betahat2 = 0.1496868;
% parms2 = [7830.103 0.007306493 0.5691455];
% loglik2 = loglikfast(zarray,observedmap,i1,i2array,D,betahat2,parms2,nu,J,0,0);
% t2 = 1.23;
% llINLA2 = 70800.23;
% llINLA2 - loglik2
% fprintf('%6.1f \n',loglik0-loglik2)
% fprintf('%4.4f & %6.2f & %7.4f & %7.4f \n',betahat2,sqrt(parms2(1)),sqrt(parms2(2)),sqrt(parms2(3)/parms2(1)))
% 
% % max.edge = 4
% betahat3 = 0.1490246;
% parms3 = [6244.351 0.01267662 0.3657674];
% loglik3 = loglikfast(zarray,observedmap,i1,i2array,D,betahat3,parms3,nu,J,0,0);
% t3 = 3.23;
% llINLA3 =  71829.78;
% llINLA3 - loglik3
% fprintf('%6.1f \n',loglik0-loglik3)
% fprintf('%4.4f & %6.2f & %7.4f & %7.4f \n',betahat3,sqrt(parms3(1)),sqrt(parms3(2)),sqrt(parms3(3)/parms3(1)))
% 
% % max.edge = 2
% betahat4 = 0.147112;
% parms4 = [3385.925 0.02939593 0.1243893];
% loglik4 = loglikfast(zarray,observedmap,i1,i2array,D,betahat4,parms4,nu,J,0,0);
% t4 = 33.2;
% llINLA4 = 72962.30;
% llINLA4 - loglik4
% fprintf('%6.1f \n',loglik0-loglik4)
% fprintf('%4.4f & %6.2f & %7.4f & %7.4f \n',betahat4,sqrt(parms4(1)),sqrt(parms4(2)),sqrt(parms4(3)/parms4(1)))
% 
% 
% 
% 
% % fit using a block-independent likelihood
% datamattemp = dlmread('readseadatavec.txt',' ',1,0);
% datamat = datamattemp(:,[4 7 10]);
% coords = datamat(:,1:2);
% z = datamat(:,3);
% 
% 
% lpoints = [1 50; 420 300];
% slp1 = diff(lpoints(:,1))/diff(lpoints(:,2))
% th = atan(1/slp1);
% rotmat = [cos(th),-sin(th);sin(th),cos(th)];
% coordstrans = ( rotmat*coords' )';
% 
% cutpoints = [0;110;145;180;210;240;270;300;325;350;375;397;422;450;520];
% plot(coordstrans(:,1),coordstrans(:,2),'.')
% hold on
% for j = 1:length(cutpoints)
%     plot([-100 200],[cutpoints(j) cutpoints(j)],'color','k')
% end
% hold off
% 
% blocks = cell(length(cutpoints)-1,1);
% npoints = zeros(length(cutpoints)-1,1);
% allinds = (1:size(coords,1))';
% for j = 2:length(cutpoints)
%    inds = coordstrans(:,2) <= cutpoints(j) & coordstrans(:,2) > cutpoints(j-1);
%    blocks{j-1} =  allinds(inds);
%    npoints(j-1) = sum( coordstrans(:,2) <= cutpoints(j) & coordstrans(:,2) > cutpoints(j-1) );
% end
% npoints
% 
% nu = 1;
% J = 3;
% 
% 
% options = optimset('display','iter','tolFun',1e-4);
% D0 = ones(length(z),1);
% funtomax = @(x) -loglikblocks(z,coords,blocks,D0,x(3),[1 exp(x(1:2))],nu,J);
% tic
% [maxparmsb,~] = fminunc(funtomax,[log(parms3(2:3)) .1468],options);
% [~,tausq] = loglikblocks(z,coords,blocks,D0,0.25,[1 exp(maxparmsb)],nu,J);
% toc
% loglikb = loglikfast(zarray,observedmap,i1,i2array,D,0.15,[tausq exp(maxparmsb(1:2))],nu,J,0,0);
% fprintf('%6.1f \n',loglik0-loglikb)
% sqrt([tausq, exp(maxparmsb(1)), exp(maxparmsb(2))/tausq, maxparmsb(3)^2])
% 
% 
% 
% % use vecchia likelihood
% cutpoints = [0,80:15:185,198:12:294,304:10:424,436,452,474, 520]';
% blocks = cell(length(cutpoints)-1,1);
% npoints = zeros(length(cutpoints)-1,1);
% allinds = (1:size(coords,1))';
% for j = 2:length(cutpoints)
%    inds = coordstrans(:,2) <= cutpoints(j) & coordstrans(:,2) > cutpoints(j-1);
%    blocks{j-1} =  allinds(inds);
% end
% nu = 1;
% J = 3;
% 
% options = optimset('display','iter','tolFun',1e-4);
% D0 = ones(length(z),1);
% funtomax = @(x) -loglikvecchia(z,coords,blocks,D0,x(3),[1 exp(x(1:2))],nu,J);
% tic
% [maxparmsb,~] = fminunc(funtomax,[log(parms3(2:3)) 0.15],options);
% [~,tausq] = loglikvecchia(z,coords,blocks,D0,0.25,[1 exp(maxparmsb)],nu,J);
% toc
% loglikb = loglikfast(zarray,observedmap,i1,i2array,D,0.15,[tausq exp(maxparmsb(1:2))],nu,J,0,0);
% fprintf('%6.1f \n',loglik0-loglikb)
% sqrt([tausq exp(maxparmsb(1)) exp(maxparmsb(2))/tausq maxparmsb(3)^2])
% 
% 



