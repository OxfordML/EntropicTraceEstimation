function [predicted_labels,D2,time,HH_center] = nn(X_train, y_train, X_test, opt)

% Nearest Neighbor classification

tStart = tic;

% PCA
if opt.pca == true
    [U, s] = myPCA(X_train);
    cs = cumsum(s) / sum(s);
    nPCA = nnz(cs < opt.pcaThres) + 1;
    P = U(:,1:nPCA).';
    X_train = myPCA_apply(P, X_train);
end

% tic;
HH_train = getHH(X_train, opt);
% toc
% tic
% HH_train2 = mex_getHH(X_train, opt.H_rows);
% toc

unique_classes = unique(y_train);
n_classes = length(unique_classes);

if strcmp(opt.metric,'SubspaceAngleFast')
    for i = 1:length(HH_train) 
        [U1,S1,V1] = svd(HH_train{i});
        s1 = diag(S1);
        ind1 = cumsum(s1)/sum(s1) < opt.SA_thr;
        ind1 = [true;ind1(1:end-1)];
        HH_train{i} = U1(:,ind1);
    end
end

if strcmp(opt.metric,'binlong') || strcmp(opt.metric,'SubspaceAngle') ||...
        strcmp(opt.metric,'SubspaceAngleFast')
    D = HHdist(HH_train, HH_train, opt);
    centerInd = findCenters(D, y_train);
    HH_center = HH_train(centerInd);
elseif strcmp(opt.metric,'JBLD') || strcmp(opt.metric,'JLD_denoise')
    HH_center = cell(1, n_classes);
    %         cparams(1:n_classes) = struct ('prior',0,'alpha',0,'theta',0);
    for ai = 1:n_classes
        X_tmp = HH_train(y_train==unique_classes(ai));
        HH_center{ai} = steinMean(cat(3,X_tmp{1:end}));
%                     HH_center{ai} = incSteinMean(cat(3,X_tmp{1:end}));
        %             d = HHdist(HH_center(ai),X_tmp,opt.metric);
        %             d(abs(d)<1e-6) = 1e-6;
        % %             phat = gamfit(d);
        %             phat = mle(d,'pdf',@gampdf,'start',[1 1],'lowerbound',[0 0],'upperbound',[1.5 inf]);
        %             cparams(ai).alpha = min(100,phat(1));
        %             if isinf(cparams(ai).alpha), keyboard;end
        %             cparams(ai).theta = max(0.01,phat(2));
        %             cparams(ai).prior = length(X_tmp) / length(X_train);
        fprintf('processed %d/%d\n',ai,n_classes);
    end
elseif strcmp(opt.metric,'AIRM')
    HH_center = cell(1, n_classes);
    for ai = 1:n_classes
        X_tmp = HH_train(y_train==unique_classes(ai));
        HH_center{ai} = karcher(X_tmp{1:end});
    end
elseif strcmp(opt.metric,'LERM')
    HH_center = cell(1, n_classes);
    for ai = 1:n_classes
        X_tmp = HH_train(y_train==unique_classes(ai));
        HH_center{ai} = logEucMean(X_tmp{1:end});
    end
elseif strcmp(opt.metric,'KLDM')
    HH_center = cell(1, n_classes);
    for ai = 1:n_classes
        X_tmp = HH_train(y_train==unique_classes(ai));
        HH_center{ai} = jefferyMean(X_tmp{1:end});
    end
end

time.trainTime = toc(tStart);

% test NN
tStart = tic;

if opt.pca == true
    X_test = myPCA_apply(P, X_test);
end

HH_test = getHH(X_test, opt);

if strcmp(opt.metric,'SubspaceAngleFast')
    for i = 1:length(HH_test) 
        [U1,S1,V1] = svd(HH_test{i});
        s1 = diag(S1);
        ind1 = cumsum(s1)/sum(s1) < opt.SA_thr;
        ind1 = [true;ind1(1:end-1)];
        HH_test{i} = U1(:,ind1);
    end
end

D2 = HHdist(HH_center,HH_test,opt);
[~,ind] = min(D2);
predicted_labels = unique_classes(ind);

time.testTime = toc(tStart);

%         % test gamma voting
%         D2 = HHdist(HH_center,X_test,opt.metric);
%         P2 = zeros(size(D2));
%         for ai = 1:size(D2,1)
%             P2(ai,:) = gampdf(D2(ai,:),...
%                 cparams(ai).alpha, cparams(ai).theta);
%         end
%         [~,ind] = max(P2);
%         predicted_labels = unique_classes(ind);



end