% Timing experiment to compare approximate and
% exact likelihoods

% set up the defaults for the plots
set(0,'defaultaxesfontsize',14)
set(0,'defaulttextfontsize',14)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'defaulttextfontname','Times New Roman')


%cd \\stat.ad.ncsu.edu\Redirect\guinness\Documents\research\carlik\code\codeforjcgs


% grid sizes and nu parameters
nvec = 100:100:300;
nuvec = [0 1];

% amount of oversampling to compute covariances
J = 2;

% parameters to use
tausq = 1;
kappasq = 1/10;
nug = 0.0;
parms = [tausq kappasq nug];

% array to save timing results
timemat = zeros(length(nvec),3,length(nuvec));

% loop over values of nu
for t1 = 1:length(nuvec);
    nu = nuvec(t1);
    
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
    for nn = 1:length(nvec)
        nn
        n1 = nvec(nn);
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
        
        
        time1 = clock;
        lq(nn, t1) = loglikQ(zarray,observedmap,i11,i21array,D,beta,parms(1:2),nu,getQ);
        time2 = clock;
        lq_(nn, t1) = loglikQ_(zarray,observedmap,i11,i21array,D,beta,parms(1:2),nu,getQ);
        %loglikfast(zarray,observedmap,i13,i23array,D,beta,parms,nu,J,0,0);
        time3 = clock;
        %loglikfast(zarray,observedmap,i13,i23array,D,beta,[parms(1:2) 0.01],nu,J,0,0);
        time4 = clock;
         
        timemat(nn,:,t1) = [etime(time2,time1),etime(time3,time2),etime(time4,time3)];
        
    end
end

% print out the results
fprintf('$ %3d^2 $ & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\ \n',[nvec' timemat(:,:,1) timemat(:,:,2)]') 

