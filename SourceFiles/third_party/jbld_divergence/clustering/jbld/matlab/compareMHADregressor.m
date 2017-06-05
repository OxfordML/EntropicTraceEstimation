function compareMHADregressor

addpath(genpath('../3rdParty'));
addpath(genpath('.'));
load ../skeleton_data/MHAD/MHAD_data_whole;

t1 = data{1+5*5}; % jump by subject 1
% show_skel_MHAD(t1);
t2 = data{1+5*5+55}; % jump by subject 2
% show_skel_MHAD(t2);
% plot([t1 t2]');
t3 = data{1+5*4+55*2}; % clap by subject 1
% show_skel_MHAD(t3);

nc = 3;
Ht1 = blockHankel(t1, [size(t1,1)*(size(t1,2)-nc), nc+1]);
% HSTLN
% [t_hat, eta, r] = hstln_mo(t1,nc);
% [U,S,V] = svd(Ht);
% r = V(:,end);
r1 = Ht1(:,1:end-1) \ Ht1(:,end);
% Ht = blockHankel(t, [size(t,1)*(size(t,2)-nc+1), nc]);
eta1 = Ht1 * [r1; -1];
plot(eta1);
norm(eta1) / length(eta1)


% t2 = data{1+5};
Ht2 = blockHankel(t2, [size(t2,1)*(size(t2,2)-nc), nc+1]);
eta2 = Ht2 * [r1; -1];
plot(eta2);
norm(eta2) / length(eta2)


Ht3 = blockHankel(t3, [size(t3,1)*(size(t3,2)-nc), nc+1]);
eta3 = Ht3 * [r1; -1];
plot(eta3);
norm(eta3) / length(eta3)


end