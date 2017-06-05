n = 5;
rng(0);
R1 = rand(n);
A = R1 * R1';
R2 = rand(n);
B = R2 * R2 + A;

JBLD(A, B)