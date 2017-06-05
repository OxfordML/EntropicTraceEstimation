function XX = getXX(X)

n = length(X);
XX = cell(1, 1, n);
for i = 1:n
    XX{i} = X{i} * X{i}.';
end

end