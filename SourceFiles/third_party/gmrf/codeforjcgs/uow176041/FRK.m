%%%%%%%%%%%%%%%   fixed rank kriging   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pred sig2FRK]=FRK(data,pred_locs,K,sigxi,sige,V,BStype,V2)

if nargin<8, V2=sparse(1:length(data),1:length(data),1); end

% observed data
lon=data(:,2);
lat=data(:,1);
z=data(:,3);

% coordinates of prediction locations
lon_pred=pred_locs(:,2); 
lat_pred=pred_locs(:,1);

n=size(z,1);  % number of observations
m=length(lon_pred); % number of prediction locations


% load basis function centers
level1 = dlmread('level1.csv', ',','A1..B32');
level2 = dlmread('level2.csv', ',','A1..B92');
level3 = dlmread('level3.csv', ',','A1..B272');

% BFs evaluated at observed locations
So=Create_S(data(:,[2 1]),level1,level2,level3,BStype);

% BFs evaluated at prediction locations
Sp=Create_S([lon_pred lat_pred],level1,level2,level3,BStype); 



%%%%%%%%%%  smooth the data  %%%%%%%%%%%%%%%

% diagonal error variance matrix
D=sigxi*V2 + sige*V;
DInv=inv(D);

% calculate r x r part of the inverse of Sigma
temp=inv(inv(K)+So'*DInv*So);

% indicator matrix for observed locations
E=sparse(m,n);
for i=1:n,
   E((lon_pred==lon(i) & lat_pred==lat(i)),i)=1;
end;

% predict smoothed residuals
SigInvZ=DInv*z-DInv*So*(temp*So'*DInv*z);
pred=Sp*(K*So'*SigInvZ)+sigxi*E*SigInvZ;



%%%%%%%%  calculate the FRK variance  %%%%%%%%

% find rows that correspond to observed locations
index=zeros(m,1);
[indrow indcol]=find(E==1);
for j=1:length(indrow), index(indrow(j),1)=indcol(j); end

% part 1
SpK=Sp*K;
p1=zeros(m,1);
for i=1:m,
    p1(i,1)=SpK(i,:)*Sp(i,:)';
end

% part 3
KSSigInv=K*So'*DInv-(K*So'*DInv*So*temp)*So'*DInv;
SpKSSigInvSK=Sp*(KSSigInv*So*K);
SigInv1=DInv*So*temp;
SigInv2=So'*DInv;
p3=zeros(m,1);
for i=1:m,
    p3(i,1)=SpKSSigInvSK(i,:)*Sp(i,:)'; 
    ind=index(i,1);
    if ind>0,
        ESigInvE=DInv(ind,ind)-SigInv1(ind,:)*SigInv2(:,ind);
        p3(i,1)=p3(i,1)+2*sigxi*Sp(i,:)*KSSigInv(:,ind)+sigxi^2*ESigInvE; 
    end;
end


%%% put the pieces of variance together
sig2FRK=p1+sigxi-p3;
