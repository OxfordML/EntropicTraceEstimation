function [logdet_me] = maxent_logdet(A)

    %%MAXENT
    max_e = max(sum(abs(A),2));
    S = A./max_e;
    d = 100;
    z = randn(size(A,1),d);
    l = 5;
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
    minA = max([1e-6, min(diag(S)*2-sum(abs(S), 2))]);
    if minA ==1
        logdet_me = log(max_e)*size(A,1);
        return
    end
    steps = (1 - 1e-2)/100;
    x = minA:steps:1;
    [p, ~, ~] = stableMaxEnt(E, x);
    r = (sum(p.*log(x')))/sum(p);
    logdet_me = log(max_e)*size(A,1) + r*size(A,1);
    %%END
end