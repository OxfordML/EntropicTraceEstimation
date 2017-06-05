% test find translation

load ../skeleton_data/MHAD/MHAD_data_whole;

t = data{1}(1:6,20:140);
plot(t');
nr = 10;
nc = size(t, 2)-nr+1;
epsilon = 1;
maxIter = 30;
iter = 1;
sigma0 = inf;

Ht = blockHankel(t, [nr nc]);
s = svd(Ht);

cvx_begin
cvx_solver mosek
variables n(size(t));
variable m(size(t,1),1);
variables a(size(t));
N = blockHankel(n, [nr nc]);
mv = kron(m, ones(1, size(t,2)));
M = blockHankel(mv, [nr nc]);
Y = Ht - N - M;
Y * c == 0; % nonconvex, it's better to try atom instead
cvx_end;

N = blockHankel(n, [nr size(t, 2)-nr+1]);
mv = kron(m, ones(1, size(t,2)));
t_clean = t - n - mv;
T = blockHankel(t_clean, [nr size(t, 2)-nr+1]);
plot([t;t-n]);
