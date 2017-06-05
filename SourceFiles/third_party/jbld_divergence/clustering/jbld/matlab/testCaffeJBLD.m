% test caffe + JBLD metric on UCF

dbstop if error

addpath(genpath('3rdParty'));

clear;close all;

dataPath = fullfile('~','research','data','UCF101','CaffeFeat');
data = parseUCF101Caffe(dataPath);

% PCA
% allData = cat(2, data{:});
% [U,mu] = pca(allData);
k = 20;
pcaData = cell(1, length(data));
for i = 1:length(data)
    [U,mu] = pca(data{i});
    pcaData{i} = pcaApply( data{i}, U, mu, k );
end

opt.H_structure = 'HHt';
% opt.H_structure = 'HtH';
% opt.metric = 'binlong';
opt.metric = 'JBLD';
% opt.metric = 'AIRM';
% opt.metric = 'LERM';
% opt.metric = 'KLDM';
% opt.Hsize = 3;
% opt.H_rows = 30;
opt.H_rows = 2;
opt.sigma = 1e-4;

v = getVelocity(pcaData);
G = getHH(v, opt);
D = HHdist(G,[],opt);
imagesc(D); colorbar;

fileName = 'D_binlong_UCF.mat';
load(fullfile(dataPath,fileName));

load (fullfile(dataPath, 'gt_labels.mat'));

% index = find(gt_labels==72 | gt_labels==95);
% D2 = D;
% D2(index, :) = [];
% D2(:, index) = [];
% gt_labels2 = gt_labels;
% gt_labels2(index) = [];

nCluster = 101;
rng(0);
[label,W] = ncutD(D,nCluster,100);

newLabel = zeros(size(label));
uniLabel = unique(label);
for i = 1:length(uniLabel)
    ind = find(label == uniLabel(i));
    newLabel(ind) = mode(gt_labels(ind));
end

accuracy = nnz(newLabel==gt_labels) / length(gt_labels);
accuracy

C = confusionMatrix(gt_labels, newLabel);
imagesc(C); colorbar;
c = diag(C);
find(c < 0.8)