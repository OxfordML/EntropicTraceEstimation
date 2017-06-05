% Timing experiment to compare approximate and
% exact likelihoods

% set up the defaults for the plots
set(0,'defaultaxesfontsize',14)
set(0,'defaulttextfontsize',14)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'defaulttextfontname','Times New Roman')


%cd /Users/jsguinne/Documents/research/carlik/code/ver5
%cd \\stat.ad.ncsu.edu\Redirect\guinness\Documents\research\carlik\code\ver5

nsims = 1;
diffvec = zeros(nsims,2);
for j = 1:nsims,
% grid sizes and nu parameters
nvec = 100;
nuvec = [0];

% amount of oversampling to compute covariances
J = 4;

% parameters to use
tausq = 1;
kappasq = 1/10;
nug = 0.01;
parms = [tausq kappasq nug];


    nu = nuvec(1);
    
    % we compute three likelihoods for each (parameter,gridsize)
    % combination--one with the approximation and no adjustment
    % to the precision matrix, one with the approximation 
    % and a periodic adjustment to the precision matrix
    % and one exact likelihood.
    getQ = @getQ_none;
    pattern1 = @getSparsePattern;
    pattern2 = @getSparsePattern_periodic;
    pattern3 = @getSparsePattern;
    
    
    % loop over grid sizes
        
        n1 = nvec(1);
        n2 = n1;
        n = n1*n2;
        boldn = [n1,n2];
        observedmap = zeros(n1,n2);
        observedmap(:) = 1;
        zarray = simulateMRF(boldn,parms,nu,2*J);
        D = ones(n,1);
        beta = 0;
        
        [i11,i21array] = pattern1(observedmap,nu);
        [i12,i22array] = pattern2(observedmap,nu);
        [i13,i23array] = pattern3(observedmap,nu);
        i1 = i13;
        i2array = i23array;
        
      
%        getQ_exact(observedmap,i1,i2array,parms,nu,zarray);
        
%         time1 = clock;
%         loglikQ(zarray,observedmap,i11,i21array,D,beta,parms(1:2),nu,getQ);
%         time2 = clock;
%         %loglikQ(zarray,observedmap,i12,i22array,D,beta,parms(1:2),nu,getQ);
%         time3 = clock;
%         loglikfast(zarray,observedmap,i13,i23array,D,beta,parms,nu,J,0,0);
%         time4 = clock;
         [L1,D1,F1] = loglikK(zarray,observedmap,parms,nu,J);
         [L2,~,~,D2,F2] = loglikfast(zarray,observedmap,i13,i23array,D,beta,parms,nu,J,0,0);
         L3 = loglikQ(zarray,observedmap,i13,i23array,D,beta,parms,nu,J);

diffvec(j,1) = L2 - L1;
diffvec(j,2) = L3 - L1;


end
diffvec

      %         time5 = clock;
%         
%         
%         timemat(nn,:,t1) = [etime(time2,time1),etime(time4,time3),etime(time5,time4)];
%         
% 
% 
% 
% pc1 = '-+';
% pc2 = '-o';
% col1 = 0.5*[1 1 1];
% col2 = 'k';
% 
% % plot the results
% for t1=1:2
% subplot(1,1,1)    
% plot(nvec(:)',(timemat(:,1,t1)).^(2/2),pc1,'color',col1)
% hold on
% plot(nvec(:)',(timemat(:,2,t1)).^(2/2),pc2,'color',col2)
% plot(nvec(:)',(timemat(:,3,t1)).^(2/2),pc2,'color',col2)
% end
% hold off
% 
% % print out the results
% timemat
% fprintf('$ %3d^2 $ & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\ \n',[nvec' timemat(:,:,1) timemat(:,:,2)]') 
% 
