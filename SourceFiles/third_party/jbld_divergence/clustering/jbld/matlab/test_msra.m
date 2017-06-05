%% test MHAD

clear;clc;close all;

addpath(genpath('../matlab'));
addpath(genpath('../3rdParty'));

load ../expData/MSRAction3D_data_AS3.mat

t1 = data{1}(1:3,:);
t2 = data{6}(1:3,:);

nc1 = 10;
nr1 = 3*(size(t1,2)-nc1+1);
Ht1 = hankel_mo(t1,[nr1 nc1]);
HHt1 = Ht1' * Ht1;
HHt1 = HHt1 / norm(HHt1,'fro');
nc2 = 10;
nr2 = 3*(size(t2,2)-nc2+1);
Ht2 = hankel_mo(t2,[nr2 nc2]);
HHt2 = Ht2' * Ht2;
HHt2 = HHt2 / norm(HHt2,'fro');

%%
% nc = 5;
% k=10; opt.metric='JLD';
% hFeat = zeros(20*k,length(data));
% for j=1:20
%     fprintf('Processing Joint %d ...\n',j)
%     HH = cell(1,length(data));
%     for i=1:length(data)
%         t = diff(data{i}(3*(j-1)+1:3*j,:),[],2);
%         nr = 3*(size(t,2)-nc+1);
%         if nr<=1,keyboard;end
%         Ht = hankel_mo(t,[nr nc]);
%         HHt = Ht' * Ht;
%         HHt = HHt / norm(HHt,'fro');
%         %     HHt = t * t';
%         I = 1e-6*eye(size(HHt));
%         HH{i} = HHt + I;
%     end
%     [label,HH_center,D,param] = ncutJLD(HH(idx_train),k,opt);
%     D2 = HHdist(HH_center,HH,opt.metric);
%     for i=1:k
%         hFeat(k*(j-1)+i,:) = gampdf(D2(i,:),param(i).alpha,param(i).theta);
%     end
% end
%%
Hsize = 540;
opt.metric='JLD';opt.H_structure = 'HHt';
HH = cell(1,length(data));
for i=1:length(data)
    t = diff(data{i},[],2);
    if strcmp(opt.H_structure,'HtH')
        nc = Hsize;
        nr = size(t,1)*(size(t,2)-nc+1);
        if nr<1, error('hankel size is too large.\n'); end
        Ht = hankel_mo(t,[nr nc]);
        HHt = Ht' * Ht;
    elseif strcmp(opt.H_structure,'HHt')
        nr = floor(Hsize/size(t,1))*size(t,1);
        nc = size(t,2)-floor(nr/size(t,1))+1;
        if nc<1, error('hankel size is too large.\n'); end
        Ht = hankel_mo(t,[nr nc]);
        HHt = Ht * Ht';
    end
    HHt = HHt / norm(HHt,'fro');
    %     HHt = t * t';
    I = 0.9*eye(size(HHt));
    HH{i} = HHt + I;
end
% get centers
X_train = HH(idx_train);
y_train = label_act(idx_train);
X_test = HH(idx_test);
y_test = label_act(idx_test);
unique_classes = unique(label_act);
n_classes = length(unique_classes);
HH_center = cell(1, n_classes);
for i = 1:n_classes
    X_tmp = X_train(y_train==unique_classes(i));
    HH_center{i} = karcher(X_tmp{1:end});
    fprintf('processed %d/%d\n',i,n_classes);
end
% test
D2 = HHdist(HH_center,X_test,opt.metric);
[~,ind] = min(D2);
predicted_labels = unique_classes(ind);
accuracy = nnz(y_test==predicted_labels)/ length(y_test);
accuracy

%%
% scale data
mx = max(hFeat,[],2); mn = min(hFeat,[],2);
hFeat = bsxfun(@rdivide,bsxfun(@minus,hFeat,(0.5*mx+0.5*mn)),0.5*mx-0.5*mn);

X_train = hFeat(:,idx_train);
y_train = label_act(idx_train);
X_test = hFeat(:,idx_test);
y_test = label_act(idx_test);

% Cind = -10:10;
% G = 10.^Cind;
% C = 2.^Cind;
% accuracyMat = zeros(length(G),length(C));
% for gi = 1:length(G)
% for ci = 1:length(C)
ly = unique(y_train);
svmModel = cell(1,length(ly));
accuracy = zeros(1,length(ly));
total_prob = zeros(length(ly), length(y_test));
for i=1:length(ly)
    y_train2 = y_train;
    y_train2(y_train==ly(i)) = 1;
    y_train2(y_train~=ly(i)) = -1;
%     y_test2 = y_test;
%     y_test2(y_test==ly(i)) = 1;
%     y_test2(y_test~=ly(i)) = -1;
    model = svmtrain(y_train2',X_train',sprintf('-s 0 -t 0 -c %d -w1 10 -w-1 1 -b 1',10));
    [predict_label, ~, prob] = svmpredict(y_test', X_test', model);
%     accuracy = nnz(predict_label==y_test2')/length(y_test2);
    %         errInd{i} = find(predict_label~=y_test2');
    %         [predict_label, ~, prob_estimates] = svmpredict(y_train2', X2_train', model);
%     accuracy(i) = nnz(predict_label==y_test2')/length(y_test2);
    total_prob(i,:) = prob';
    svmModel{i} = model;
end

[~,ind] = max(total_prob);
final_predicted_labels = ly(ind);

total_accuracy = length(find(y_test == final_predicted_labels))...
    / length(y_test);

%     accuracy
    fprintf('\naccuracy is %f\n',total_accuracy);
%     accuracyMat(gi,ci) = mean(accuracy);
% end
% end

%%
D = HHdist(HH,HH,'JLD');
imagesc(D(1:55,1:55))

