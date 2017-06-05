% plot the example lattice

% set up the defaults for the plots
set(0,'defaultaxesfontsize',14)
set(0,'defaulttextfontsize',14)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'defaulttextfontname','Times New Roman')


cd \\stat.ad.ncsu.edu\Redirect\guinness\Documents\research\carlik\code\codeforjcgs





n1 = 24;
n2 = 12;

observedinds = ones(n1,n2);
observedinds(7,5) = 0;
observedinds(8,5) = 0;
observedinds(18,8) = 0;

neighborinds = ones(n1,n2);
neighborinds(1,:) = 0;
neighborinds(:,1) = 0;
neighborinds(n1,:) = 0;
neighborinds(:,n2) = 0;

gr = 0.7*[1 1 1];

figure('position',[400 400 600 300])
subplot('position',[0.01 0.02 0.95 0.9])
for j1 = 1:n1,
    for j2 = 1:n2,
        if observedinds(j1,j2) == 1
            plot( j1, j2 ,'o','color','k','markersize',10)
            hold on
        else
           neighborinds(j1,j2) = 0;
           neighborinds(j1+1,j2) = 0;
           neighborinds(j1-1,j2) = 0;
           neighborinds(j1,j2+1) = 0;
           neighborinds(j1,j2-1) = 0;
        end
    end
end
for j1 = 1:n1,
    for j2 = 1:n2,
        if neighborinds(j1,j2) == 1
            plot( j1, j2 ,'.','color',gr,'markersize',20)
            hold on
        end
    end
end
hold off
axis off
