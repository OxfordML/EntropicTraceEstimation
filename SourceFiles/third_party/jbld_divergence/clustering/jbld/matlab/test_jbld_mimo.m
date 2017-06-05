
clear all
clc
close all
ny = 2; % dimension of the output
Lsys = 2; % number of systems
nx = 4; % dimension of the states
As = cell(Lsys,1);
Cs = cell(Lsys,1);
for i = 1:Lsys
%     As{i} = RandOrthMat(nx);
    if i==1, As{i} = RandOrthMat(nx);
    else As{i} = As{1}; end
%      Cs{i} = randn(ny,size(As{i},1));
    Cs{i} = zeros(ny,size(As{i},1));
    Cs{i}(:,1:ny) = eye(ny);
end


Lseq = 10;

% C = zeros(ny,nx);C(:,1:ny) = eye(ny);

Dens = [];
for sysId = 1:Lsys
    %Bs{sysId} = zeros(size(As{sysId},1),1); Bs{sysId}(1) = 1;
    Bs{sysId} = rand(size(As{sysId},1),1);
    nu = size(Bs{sysId},2);

    sys = ss(As{sysId},Bs{sysId},Cs{sysId},zeros(ny,nu),-1);
    [~,den] = tfdata(tf(sys),'v');
    Dens{sysId} = den{1};
end

%%
L = 30;
Xs = cell(Lsys,Lseq);
Ys = cell(Lsys,Lseq);
for sysId = 1:Lsys   
    for seqId = 1:Lseq
        nx = size(As{sysId},1);
        Xs{sysId,seqId} = zeros(nx,L);
        Xs{sysId,seqId}(:,1) = randn(nx,1);
        Ys{sysId,seqId} = zeros(ny,L);
        for i = 1:L
            Xs{sysId,seqId}(:,i+1) = As{sysId}*Xs{sysId,seqId}(:,i);
            Ys{sysId,seqId}(:,i) = Cs{sysId}*Xs{sysId,seqId}(:,i);
        end
    end
end
%save('data_jbld_mimo');
%}
%%

Hs = cell(Lsys,Lseq);
HHts = cell(Lsys,Lseq);
HtHs = cell(Lsys,Lseq);
nHHts = cell(Lsys,Lseq);
nHtHs = cell(Lsys,Lseq);

colH = size(As{1},1)+1;
for sysId = 1:Lsys
    for seqId = 1:Lseq
%         Hs{sysId,seqId} = blkhankel(Ys{sysId,seqId},colH);
        Hs{sysId,seqId} = blockHankel(Ys{sysId,seqId},[colH, size(Ys{sysId,seqId},2)-colH+1]);
        HHts{sysId,seqId} = Hs{sysId,seqId}*transpose(Hs{sysId,seqId});
        HtHs{sysId,seqId} = transpose(Hs{sysId,seqId})*Hs{sysId,seqId};
        
        nHHts{sysId,seqId} = HHts{sysId,seqId}/norm(HHts{sysId,seqId},'fro');
        nHtHs{sysId,seqId} = HtHs{sysId,seqId}/norm(HtHs{sysId,seqId},'fro');
    end
end
%% single values
ssHHt = cell(1,Lsys);
ssHtH = cell(1,Lsys);

for sysId = 1:Lsys
    ssHHt{sysId} = [];
    ssHtH{sysId} = [];
    for seqId = 1:Lseq
        ssHHt{sysId} = [ssHHt{sysId},svd(HHts{sysId,seqId})];
        ssHtH{sysId} = [ssHtH{sysId},svd(HtHs{sysId,seqId})];
    end
end
%{
resx = cell(Lsys,1);
resy = cell(Lsys,1);
for denId = 1:Lsys
    resx{denId} = zeros(Lseq,Lsys);
for sysId = 1:Lsys
    for seqId = 1:Lseq
        xx = Dens{denId}*Hs{sysId,seqId}(1:2:end,:);
        resx{denId}(seqId,sysId) = xx*xx';
        yy = Dens{denId}*Hs{sysId,seqId}(2:2:end,:);
        resy{denId}(seqId,sysId) = yy*yy';
    end
end
end
%}
%% distance 
dist_M_HHt = zeros(Lsys*Lseq);
dist_M_HtH = zeros(Lsys*Lseq);
dist_M_nHHt = zeros(Lsys*Lseq);
dist_M_nHtH = zeros(Lsys*Lseq);

dist_M_nHHt_bl = zeros(Lsys*Lseq);
dist_M_nHtH_bl = zeros(Lsys*Lseq);

nHHts = nHHts';
nHtHs = nHtHs';
HHts = HHts';
HtHs = HtHs';

regul = 1e-4;
for i = 1:Lsys*Lseq
    for j = 1:Lsys*Lseq
%         dist_M_nHHt(i,j) = dist_JBLD(nHHts{i}+regul*eye(size(nHHts{i})),nHHts{j}+regul*eye(size(nHHts{j})));
%         dist_M_nHtH(i,j) = dist_JBLD(nHtHs{i}+regul*eye(size(nHtHs{i})),nHtHs{j}+regul*eye(size(nHtHs{j})));
%         
%         dist_M_HHt(i,j) = dist_JBLD(HHts{i}+regul*eye(size(HHts{i})),HHts{j}+regul*eye(size(HHts{j})));
%         dist_M_HtH(i,j) = dist_JBLD(HtHs{i}+regul*eye(size(HtHs{i})),HtHs{j}+regul*eye(size(HtHs{j})));
        
        dist_M_nHHt(i,j) = JBLD(nHHts{i}+regul*eye(size(nHHts{i})),nHHts{j}+regul*eye(size(nHHts{j})));
        dist_M_nHtH(i,j) = JBLD(nHtHs{i}+regul*eye(size(nHtHs{i})),nHtHs{j}+regul*eye(size(nHtHs{j})));
        
        dist_M_HHt(i,j) = JBLD(HHts{i}+regul*eye(size(HHts{i})),HHts{j}+regul*eye(size(HHts{j})));
        dist_M_HtH(i,j) = JBLD(HtHs{i}+regul*eye(size(HtHs{i})),HtHs{j}+regul*eye(size(HtHs{j})));
        
        dist_M_nHHt_bl(i,j) = 2-norm(nHHts{i}+nHHts{j},'fro');
        dist_M_nHtH_bl(i,j) = 2-norm(nHtHs{i}+nHtHs{j},'fro');
    end
end

figure,subplot(2,3,1),imagesc(dist_M_HHt),colorbar,title('HHt_JBLD');
subplot(2,3,4),imagesc(dist_M_HtH),colorbar,title('HtH_JBLD');

subplot(2,3,2),imagesc(dist_M_nHHt),colorbar,title('nHHt_JBLD');
subplot(2,3,5),imagesc(dist_M_nHtH),colorbar,title('nHtH_JBLD');

subplot(2,3,3),imagesc(dist_M_nHHt_bl),colorbar,title('nHHt_BL');
subplot(2,3,6),imagesc(dist_M_nHtH_bl),colorbar,title('nHtH_BL');






