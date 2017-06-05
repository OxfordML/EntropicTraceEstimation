function [U, s] = myPCA(X)

XX = getXX(X);
A = cell2mat(XX);
XXtot = sum(A, 3);
[eigVec, eigVal] = eig(XXtot);
[s, ind] = sort(-sqrt(abs(diag(eigVal))));
s = -s;
U = eigVec(:, ind);

end