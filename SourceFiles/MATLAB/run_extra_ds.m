load(['path_to_dataset/',ds_name,'.mat'])
A = Problem.A;
S = A./max(sum(abs(A),2));
n = max(sum(abs(A),2));
T = 1.1093e6;
TS = logdet(S);

minA = max([min(diag(S)*2 - sum(abs(S),2)), 1e-8]);

x = minA:0.01:1;
d = 30;
    z = randn(size(A,1),d);
    l = 10;
    E = zeros(l,1);

    for i=1:d
       z_ = z(:,i)/sqrt(z(:,i)'*z(:,i));
       Ez = S*z_;
       for j=2:l
          E(j) = E(j) + ((z_'*Ez)/d);
          Ez = S*Ez;
       end
    end
    

    E(1) = 1;
    E(2) = trace(S)/size(A,1);
    E(3) = sum(sum(S.^2))/size(A,1);
[p, N, moments] = stableMaxEnt(E, x);
r = (sum(p.*log(x')))/sum(p);%mean(plnp)/mean(p_);
rel_error = abs(r*size(A,1) - TS  )/abs(TS)
