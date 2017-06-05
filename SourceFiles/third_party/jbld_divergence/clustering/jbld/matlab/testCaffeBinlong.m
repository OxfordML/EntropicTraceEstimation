% test caffe + binlong's metric on UCF

dbstop if error

addpath(genpath('3rdParty'));

% clear;close all;
% 
% dataPath = fullfile('~','research','data','UCF101');
% fileName = 'D_binlong_UCF.mat';
% load(fullfile(dataPath,fileName));

load (fullfile(dataPath, 'gt_labels.mat'));

index = find(gt_labels==72 | gt_labels==95);
D2 = D;
D2(index, :) = [];
D2(:, index) = [];
gt_labels2 = gt_labels;
gt_labels2(index) = [];

nCluster = 99;
rng(0);
[label,W] = ncutD(D2,nCluster,100);

newLabel = zeros(size(label));
uniLabel = unique(label);
for i = 1:length(uniLabel)
    ind = find(label == uniLabel(i));
    newLabel(ind) = mode(gt_labels2(ind));
end

accuracy = nnz(newLabel==gt_labels2) / length(gt_labels2);
accuracy

C = confusionMatrix(gt_labels2, newLabel);
imagesc(C); colorbar;
c = diag(C);
find(c < 0.8)