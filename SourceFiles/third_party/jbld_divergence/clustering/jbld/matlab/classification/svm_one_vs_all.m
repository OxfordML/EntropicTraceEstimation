function [total_accuracy, class_wise_accuracy, confusion_matrix] =...
    svm_one_vs_all(X_train, X_test,...
    y_train, y_test, C_val)

% % normalize each instance
% sumXr = sum(X_train);
% X_train = bsxfun(@rdivide, X_train, sumXr);
% sumXe = sum(X_test);
% X_test = bsxfun(@rdivide, X_test, sumXe);
% % normailize each feature
% maxXr = max(X_train,[],2);
% minXr = min(X_train,[],2);
% rangeXr = maxXr - minXr;
% tempXr = bsxfun(@minus, X_train, minXr);
% tempXr = bsxfun(@rdivide, tempXr, rangeXr);
% X_train = tempXr * 2 - 1;
% tempXe = bsxfun(@minus, X_test, minXr);
% tempXe = bsxfun(@rdivide, tempXe, rangeXr);
% X_test = tempXe * 2 - 1;

unique_classes = unique(y_train);
n_classes = length(unique_classes);
n_test_samples = length(y_test);
svmModel = cell(1,n_classes);
test_prediction_prob = zeros(n_classes, length(y_test));
for i=1:n_classes
    y_train2 = y_train;
    y_train2(y_train==unique_classes(i)) = 1;
    y_train2(y_train~=unique_classes(i)) = -1;
    class_imbalance_ratio = nnz(y_train2==-1) / nnz(y_train2==1);
    model = svmtrain(y_train2,X_train',sprintf('-s 0 -t 0 -c %d -w1 %f -b 1 -q',C_val,class_imbalance_ratio));
    [~, ~, prob] = svmpredict(y_test, X_test', model, '-b 1 -q');
    test_prediction_prob(i,:) = prob(:,1)';
    svmModel{i} = model;
end

[~, ind] = max(test_prediction_prob);
final_predicted_labels = unique_classes(ind);

class_wise_accuracy = zeros(n_classes, 1);
confusion_matrix = zeros(n_classes, n_classes);
for i = 1:n_classes
    temp = find(y_test == unique_classes(i));
    class_wise_accuracy(i) =...
        length(find(final_predicted_labels(temp) == unique_classes(i)))...
        / length(temp);
    
    confusion_matrix(i, :) = hist(final_predicted_labels(temp), unique_classes) / length(temp);
end

total_accuracy = length(find(y_test == final_predicted_labels))...
    / n_test_samples;

end