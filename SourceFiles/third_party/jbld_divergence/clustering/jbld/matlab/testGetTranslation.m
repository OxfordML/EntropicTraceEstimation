% test find translation

load ../skeleton_data/MHAD/MHAD_data_whole;

t = data{1}(1:2,1:30);
plot(t);
nr = 10;
nc = size(t, 2)-nr+1;
epsilon = 1;
maxIter = 30;
iter = 1;
sigma0 = inf;

Ht = blockHankel(t, [nr nc]);
W1 = eye(nr, nr);
W2 = eye(nc, nc);
s = svd(Ht);
sigma = sigma0;
while s(end) > 1e-3 && iter < maxIter
    cvx_begin
    cvx_solver mosek
    variables n(size(t));
    variable m(size(t,1),1);
    variables P(nr, nr) Q(nc, nc);
    N = blockHankel(n, [nr nc]);
    mv = kron(m, ones(1, size(t,2)));
    M = blockHankel(mv, [nr nc]);
    Y = Ht - N - M;
    [P Y;Y' Q] == semidefinite(nr+nc);
    obj = trace(W1*P) + trace(W2*Q);
    norm(n,inf) <= epsilon;
    minimize(obj)
    cvx_end
    sy = svd(Y);
    sigma = min([sy(2) sigma]);
    W1 = inv(P + sigma*eye(nr));
    W2 = inv(Q + sigma*eye(nc));
    Y = Ht - N - M;
    s = svd(Y);
    s'
    
    iter = iter + 1;
end
N = blockHankel(n, [nr size(t, 2)-nr+1]);
mv = kron(m, ones(1, size(t,2)));
t_clean = t - n - mv;
T = blockHankel(t_clean, [nr size(t, 2)-nr+1]);
plot([t;t-n]);
