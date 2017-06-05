% Simulation study for estimating CAR parameters
set(0,'defaultaxesfontsize',14)
set(0,'defaulttextfontsize',14)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'defaulttextfontname','Times New Roman')


%cd \\stat.ad.ncsu.edu\Redirect\guinness\Documents\research\carlik\code\codeforjcgs


% define the grid dimensions
n1 = 50;
n2 = 50;
boldn = [n1,n2];
n = n1*n2; % total number of grid points

% observedmap(i,j) = 1 if data are observed at location (i,j)
observedmap = zeros(n1,n2);
observedmap(:) = 1;

% In the simulation, we vary nu and kappasq
nuvec = [0 1];
kappasqvec = [5 10 20].^(-2);
% tausq is always 1
tausq = 2;
nug = 0.00;

% number of simulated data sets at each (nu,kappasq)
% parameter combination
nsims = 10;

% set up array to save parameters
allparms = zeros(nsims,4,3,2);
% amount of oversampling to compute covariances
J = 4;
% options for fminunc()
options = optimset('display','off');
D = ones(n,1);
beta = 0;        

warning off
% loop over nuvec
for t1 = 1:length(nuvec)

    nu = nuvec(t1);
    
    % set up options if nu = 0 or if nu >= 0
    if nu==0, % use the precision adjustment
        getQ = @getQ_precadj; 
        pattern1 = @getSparsePattern;
        pattern2 = @getSparsePattern;
    else  % use the periodic adjustment
        getQ = @getQ_none; 
        pattern1 = @getSparsePattern_periodic;
        pattern2 = @getSparsePattern;
    end
    
    % this tells us the sparsity pattern in the precision matrix Q
    [i11,i21array] = pattern1(observedmap,nu);
    [i1,i2array] = pattern2(observedmap,nu);    
    
    % loop over values of kappasq
    for k1 = 1:3
        
        kappasq = kappasqvec(k1);
        parms = [tausq kappasq nug];
        % array to save parameters
        svparms = zeros(nsims,4);
        
        % loop over the simulated datasets        
        for j = 1:nsims,

            % simulate the data
            zarray = simulateMRF(boldn,parms,nu,2*J);
            
            % loglikQ computes an approximate likelihood using Q
            funtomax = @(x) -loglikQ(zarray,observedmap,i11,i21array,D,beta,[exp(x(1:2)) 0],nu,getQ);
            [maxparms] = fminunc(funtomax,log(parms(1:2)),options);
            svparms(j,1:2) = exp(maxparms);
            
            
            funtomax = @(x) -loglikfast(zarray,observedmap,i1,i2array,D,beta,[1 exp(x) 0],nu,J,0,1);
            [maxparms,maxlik] = fminunc(funtomax,log(parms(2)),options);
            [~,~,sv2] = loglikfast(zarray,observedmap,i1,i2array,D,beta,[1 exp(maxparms) 0],nu,J,0,1);
            svparms(j,3:4) = [sv2 exp(maxparms)];

            fprintf('.')
            if mod(j,50)==0, fprintf(' %4d \n',j); end
        end
        
        fprintf('\n\n');
        
        allparms(:,:,k1,t1) = svparms;
        
    end
end
warning on



save('simulationresults','allparms');

