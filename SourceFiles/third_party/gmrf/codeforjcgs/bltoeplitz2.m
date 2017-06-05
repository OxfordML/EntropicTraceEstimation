function [ outmat ] = bltoeplitz2( inmat )
% takes in a 2D array of size n1 x n2 and creates a 
% block toeplitz with toeplitz block matrix from it

[n1,n2] = size(inmat);

outmat = zeros(n1*n2);
indsarray = reshape((1:n1*n2)',n1,n2);


for j = 1:n1,
    for k = 1:n2,
        col = j + (k-1)*n1;
        inds1 = abs( (-j+1 : -j+n1)' ) + 1;
        inds2 = abs( (-k+1 : -k+n2)' ) + 1;
        inds = indsarray(inds1,inds2);
        outmat(:,col) = inmat(inds(:));
    end
end

end

