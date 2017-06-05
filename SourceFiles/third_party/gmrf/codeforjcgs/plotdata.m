% Plotting AOT data
set(0,'defaultaxesfontsize',10)
set(0,'defaulttextfontsize',10)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'defaulttextfontname','Times New Roman')


%cd \\stat.ad.ncsu.edu\Redirect\guinness\Documents\research\carlik\code\codeforjcgs
load('redseadata.mat')


boldn = size(subdata);
zarray = double(subdata);
n1 = boldn(1);
n2 = boldn(2);
observedmap = double(observedmap);


cmap = colormap;
cmap(1,:) = [1 1 1];


figure('position',[600 100 600 250])
imagesc(zarray')
axis xy
axis off
colormap(cmap)
colorbar
text(50,250,'\leftarrow   N')
title('Aerosol Optical Thickness at 869 nm')


datamattemp = dlmread('readseadatavec.txt',' ',1,0);
datamat = datamattemp(:,[4 7 10]);
coords = datamat(:,1:2);


lpoints = [1 50; 420 300];
slp1 = diff(lpoints(:,1))/diff(lpoints(:,2));
th = atan(1/slp1);
rotmat = [cos(th),-sin(th);sin(th),cos(th)];
rotmat2 = [cos(th),sin(th);-sin(th),cos(th)];
coordstrans = ( rotmat*coords' )';

cutpoints1 = [0;110;145;180;210;240;270;300;325;350;375;397;422;450;520];
cutpoints2 = [0,80:15:185,198:12:294,304:10:424,436,452,474, 520]';
m1 = length(cutpoints1);
m2 = length(cutpoints2);
xytransL1 = [cutpoints1,-20*ones(m1,1)];
xytransR1 = [cutpoints1,100*ones(m1,1)];
xyL1 = (rotmat2 \ xytransL1')';
xyR1 = (rotmat2 \ xytransR1')';

xytransL2 = [cutpoints2,-20*ones(m2,1)];
xytransR2 = [cutpoints2,100*ones(m2,1)];
xyL2 = (rotmat2 \ xytransL2')';
xyR2 = (rotmat2 \ xytransR2')';


figure('position',[600 100 700 220])
subplot('position',[0.035 0.07 0.45 0.8])
imagesc(1:420,1:300,zarray')
axis xy
axis off
hold on
for j = 1:m1,
    plot([xyL1(j,1);xyR1(j,1)],[xyL1(j,2);xyR1(j,2)],'color','k')
end
hold off
colormap(cmap)
text(50,250,'\leftarrow   N')
title('Partition for Independent Blocks Likelihood Approximation')

subplot('position',[0.525 0.07 0.45 0.8])
imagesc(1:420,1:300,zarray')
axis xy
axis off
hold on
for j = 1:m2,
    plot([xyL2(j,1);xyR2(j,1)],[xyL2(j,2);xyR2(j,2)],'color','k')
end
hold off
colormap(cmap)
text(50,250,'\leftarrow   N')
title('Partition for Vecchia Likelihood Approximation')

