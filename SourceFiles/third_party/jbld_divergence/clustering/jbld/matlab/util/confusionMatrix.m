function C = confusionMatrix(gt, label)

ugt = unique(gt);
n = length(ugt);
C = zeros(n, n);
for i = 1:n
    for j = 1:n
    ind = find(gt == ugt(i));
    C(i,j) = nnz(label(ind) == ugt(j)) / length(ind);
%     h = hist(label(ind), 1:n);
%     C(i, :) = h / length(ind);
    end
end

end