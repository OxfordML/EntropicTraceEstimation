function featPCA = myPCA_apply(P, feat)

n = length(feat);
featPCA = cell(1, n);
for i = 1:n
    featPCA{i} = P * feat{i};
end

end