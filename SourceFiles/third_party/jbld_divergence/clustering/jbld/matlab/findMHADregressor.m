function findMHADregressor

dbstop if error

% do PCA on MHAD data, run HSTLN and find the regressor

load ../skeleton_data/MHAD/MHAD_data_whole;

t = data{1}(:,20:140);
plot(t');

% PCA
[coef, score] = pca(t');
t_pca = score.';

tm = mean(t, 2);
tc = bsxfun(@minus, t, tm);
[U,S,V] = svd(t);
t_pca = S*V';

% HSTLN
[t_hat, eta, r] = hstln_mo(t_pca(1:10,:),7);

plot(eta');


end