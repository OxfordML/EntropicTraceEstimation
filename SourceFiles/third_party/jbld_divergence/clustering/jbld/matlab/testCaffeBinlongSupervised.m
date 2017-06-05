% test caffe + binlong's metric on UCF, supervised classification

dbstop if error

addpath(genpath('3rdParty'));
addpath(genpath('.'));

% clear;close all;

dataPath = fullfile('~','research','data','UCF101');
% fileName = 'D_binlong_UCF.mat';
% load(fullfile(dataPath,fileName));

load(fullfile(dataPath, 'gt_labels.mat'));

load(fullfile(dataPath, 'trainIdx_split2.mat'));
load(fullfile(dataPath, 'testIdx_split2.mat'));
y_train = gt_labels(trainIdx);
y_test = gt_labels(testIdx);

D2 = D(testIdx, trainIdx);

[~, ind] = min(D2, [], 2);
predict_labels = y_train(ind);

accuracy = nnz(predict_labels==y_test) / length(y_test);
accuracy

C = confusionMatrix(y_test, predict_labels);
imagesc(C); colorbar;
c = diag(C);
find(c < 0.8)