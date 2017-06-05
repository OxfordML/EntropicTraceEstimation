% test JBLD on shared poles

addpath(genpath('.'));

clear; clc;

p1 = [0.9, 0.7];
p2 = [ 0.8, 0.6];
p3 = [ 0.9, 0.6];

c1 = -fliplr(poly(p1));
c2 = -fliplr(poly(p2));
c3 = -fliplr(poly(p3));

jbld = zeros(1,50);
binlong = zeros(1,50);
for j = 1:50

x1_0 = rand(1, 2);
x2_0 = rand(1, 2);
x3_0 = rand(1, 2);

y1 = zeros(1, 50);
y2 = zeros(1, 50);
y3 = zeros(1, 50);

y1(1:2) = x1_0;
y2(1:2) = x2_0;
y3(1:2) = x3_0;
for i = 3:50
    y1(i) = c1(1:end-1) * y1(i-2:i-1).';
    y2(i) = c2(1:end-1) * y2(i-2:i-1).';
    y3(i) = c3(1:end-1) * y3(i-2:i-1).';
end

nr = 5;
nc = length(y1) - nr + 1;
H1 = blockHankel(y1, [nr nc]);
H2 = blockHankel(y2, [nr nc]);
H3 = blockHankel(y3, [nr nc]);

sigma = 1e-4;
HH1 = H1 * H1';
HH2 = H2 * H2';
HH3 = H3 * H3';
HH1 = HH1 / norm(HH1, 'fro');
HH2 = HH2 / norm(HH2, 'fro');
HH3 = HH3 / norm(HH3, 'fro');

jbld12(j) = JBLD(HH1 + sigma * eye(nr, nr), HH2 +sigma * eye(nr, nr));
jbld13(j) = JBLD(HH1 + sigma * eye(nr, nr), HH3 +sigma * eye(nr, nr));
binlong12(j) = 2 - norm(HH1 + HH2, 'fro');
binlong13(j) = 2 - norm(HH1 + HH3, 'fro');

end

[jbld12; jbld13]
plot([jbld12; jbld13]')
legend jbld12 jbld13