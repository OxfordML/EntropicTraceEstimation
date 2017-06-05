% plot results of simulation study
set(0,'defaultaxesfontsize',12)
set(0,'defaulttextfontsize',12)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'defaulttextfontname','Times New Roman')


cd \\stat.ad.ncsu.edu\Redirect\guinness\Documents\research\carlik\code\codeforjcgs

load('simulationresults.mat')

% define the grid dimensions
n1 = 100;
n2 = 100;
boldn = [n1,n2];
n = n1*n2; % total number of grid points

% observedmap(i,j) = 1 if data are observed at location (i,j)
observedmap = zeros(n1,n2);
observedmap(:) = 1;
% set the values in this square to be missing
%observedmap(5:10,5:10) = 0;

% In the simulation, we vary nu and kappasq
nuvec = [0 1];
kappasqvec = [5 10 20].^(-2);
% tausq is always 1
tausq = 2;
nug = 0.00;
J = 4;

% number of simulated data sets at each (nu,kappasq)
% parameter combination
nsims = 100;




% plot example simulations from all six parameter combos
figure('position',[100 200 600 300])
for t1 = 1:2,
    for k1 = 1:length(kappasqvec),

        nu = nuvec(t1);
        kappasq = kappasqvec(k1);
        parms = [tausq kappasq nug];
        zarray = simulateMRF(boldn,parms,nu,2*J);

        subplot(2,3,k1+(t1-1)*3)
        imagesc(zarray.*observedmap)
        title(['\nu = ' int2str(nu) ', \kappa = 1/' int2str(1/sqrt(kappasq))])
    end
end




pc1 = '+';
pc2 = 'o';
lw1 = 10;
lw2 = 2;
col1 = 0.9*[1 1 1];
col2 = 'b';
col3 = 'k';
col4 = 0.5*[1 1 1];

% plot the results of the simulation
figure('position',[100 200 700 400])
xlmmat = [-0.9 -0.5;-1.4 -.6;-2.2 -0.5;...
          -0.9 -0.5;-1.4 -.6;-2.2 -0.5];
t1 = 1;
sq1 = allparms(:,1,:,t1);
sq2 = allparms(:,3,:,t1);
ry1 = [min([sq1(:);sq2(:)]) max([sq1(:);sq2(:)])];
t1 = 2;
sq1 = allparms(:,1,:,t1);
sq2 = allparms(:,3,:,t1);
ry2 = [min([sq1(:);sq2(:)]) max([sq1(:);sq2(:)])];
ylmmat = [ry1(1)-0.2*(ry1(2)-ry1(1)) ry1(2)+0.2*(ry1(2)-ry1(1));...
    ry2(1)-0.2*(ry2(2)-ry2(1)) ry2(2)+0.2*(ry2(2)-ry2(1))];
for t1 = 1:2
    for k1 = 1:length(kappasqvec)
        subplot(2,3,k1+(t1-1)*3)
        kappasq = kappasqvec(k1);
        nu = nuvec(t1);
        xlm = xlmmat(k1+(t1-1)*3,:);
        ylm = ylmmat(t1,:);
        plot(xlm,[tausq tausq],'linewidth',lw1,'color',col1)
        hold on
            plot(1/2*log10([kappasq kappasq]),ylm,'linewidth',lw1,'color',col1)
            plot(xlm,mean(allparms(:,1,k1,t1))*[1 1],'color',col3,'linewidth',lw2)
            plot(1/2*mean(log10(allparms(:,2,k1,t1)))*[1 1],ylm,'color',col3,'linewidth',lw2)
            plot(xlm,mean(allparms(:,3,k1,t1))*[1 1],'color',col4,'linewidth',lw2)
            plot(1/2*mean(log10(allparms(:,4,k1,t1)))*[1 1],ylm,'color',col4,'linewidth',lw2)
            plot(1/2*log10(allparms(:,2,k1,t1)),allparms(:,1,k1,t1),pc1,'color',col3)
            plot(1/2*log10(allparms(:,4,k1,t1)),allparms(:,3,k1,t1),pc2,'color',col4)
            xlim(xlm)
            ylim(ylm)
            if t1==2, xlabel('log_{10} \kappa'); end
            if k1==1, ylabel('\tau^2'); end
            title(['\nu = ' int2str(nu) ', \kappa = 1/' int2str(1/sqrt(kappasq))])
        hold off
    end
end


