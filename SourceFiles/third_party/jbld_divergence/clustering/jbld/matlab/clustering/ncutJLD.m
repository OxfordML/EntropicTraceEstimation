function [label,X_center,D,cparams] = ncutJLD(X,k,opt)
% ncutJLD:
% perform kmeans clustering on covariance matrices with JLD metric
% Input:
% X: an N-by-1 cell vector
% k: the number of clusters
% Output:
% label: the clustered labeling results

% N = length(X);
% D2 = zeros(N);
D = HHdist(X,[],opt);
D = D - min(D(:));
% load sD;
% D = sD;

kNN = 50;
D2 = D;
for j=1:size(D,1)
    [ignore,ind] = sort(D(:,j));
    D2(ind(kNN+1:end),j) = Inf;
    D2(ind(1:kNN),j) = D(ind(1:kNN),j) / max(D(ind(1:kNN),j)) * 0.5;
end
% D = min(D2,D2');%(B+B')/2;
D = (D2+D2')/2;
% D = max(D2,D2');

scale_sig = 1;
W = exp(-D.^2/(2*scale_sig^2));
NcutDiscrete = ncutW(W, k);
label = sortLabel_count(NcutDiscrete);

cparams(1:k) = struct ('alpha',0,'theta',0);
X_center = cell(1, k);
for j=1:k
    
    if nnz(label==j)==1
        X_center{j} = X{label==j};
        continue;
    elseif nnz(label==j)==0
        error('cluster is empty.\n');
    end
    
    if strcmp(opt.metric,'JBLD')
        X_center{j} = steinMean(cat(3,X{label==j}));
        %         d = HHdist(X_center(j),X(label==j),'JLD');
        %         d(abs(d)<1e-6) = 1e-6;
        %         param = gamfit(d);
        %         cparams(j).alpha = min(100,param(1));
        %         if isinf(cparams(j).alpha), keyboard;end
        %         cparams(j).theta = max(0.01,param(2));
    elseif strcmp(opt.metric,'AIRM')
        X_center{j} = karcher(X{label==j});
    elseif strcmp(opt.metric,'LERM')
        X_center{j} = logEucMean(X{label==j});
    elseif strcmp(opt.metric,'KLDM')
        X_center{j} = jefferyMean(X{label==j});
    elseif strcmp(opt.metric,'binlong')
        X_center{j} = findCenter(X(label==j));
    end
end

end

function center = findCenter(X)

n = length(X);
D = zeros(n,n);
for i=1:n
    for j=i+1:n
        %         D(i,j) = hankeletAngle(X{i},X{j},thr);
        D(i,j) = 2 - norm(X{i}+X{j},'fro');
    end
end
D = D + D';
d = sum(D);
[~,ind] = min(d);
center = X{ind};

end