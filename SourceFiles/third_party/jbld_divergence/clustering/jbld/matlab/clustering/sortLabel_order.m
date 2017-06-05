function label = sortLabel_order(NcutDiscrete, order)

% label the data from low order to high order
n = size(NcutDiscrete, 1);
k = size(NcutDiscrete, 2);
label = zeros(n, 1);
index = cell(k, 1);
order_index = zeros(k, 1);

for i = 1:k
    index{i} = find(NcutDiscrete(:, i));
    order_index(i) = mean(order(index{i}));
end

[~, I] = sort(order_index);

for i = 1:k
    label(index{I(i)}) = i;
end

end