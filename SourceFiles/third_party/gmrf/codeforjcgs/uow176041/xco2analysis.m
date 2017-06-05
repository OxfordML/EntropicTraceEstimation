%%%%%%%   Example Analysis of XCO2 Data    %%%%%%%%%%

clear

%%% load the data
load 'xco2data'
% the variables names for the measurements are z, lat, and lon
% the dataset also includes prediction coordinates predlat and predlon

%%% load the basis function centers, of three resolutions
level1 = dlmread('level1.csv', ',','A1..B32');
level2 = dlmread('level2.csv', ',','A1..B92');
level3 = dlmread('level3.csv', ',','A1..B272');

%%% create the basis-function matrix
So=Create_S([lon lat],level1,level2,level3,1);

%%% create other necessary quantities
n=length(z);
V_eps=sparse(1:n,1:n,1);


%%%%%%%%%    detrend the data  %%%%%%%%%

% create the matrix of predictors (intercept and longitude)
Xo=[ones(length(z),1) lon];

% ordinary least squares estimation of the trend coefficients
betahat=regress(z,Xo);

% detrending
z_tilde=z-Xo*betahat;


%%%%%%%  delete the tails to evaluate robustness
% tails=abs(z_tilde)>4*std(z_tilde); lon=lon(~tails); lat=lat(~tails);
% z_tilde(tails)=[]; n=length(z_tilde); V_eps=sparse(1:n,1:n,1);
% So=Create_S([lon lat],level1,level2,level3,1);



%%%%%%%%%    estimate parameters  %%%%%%%%%

%%%% measurement error

% true measurement error
sig2_eps = 0.25;

% variogram estimation of measurement-error and FSV variances
binendpoints=17:3:29;
[sig2_eps_variogram sig2_xi_variogram variogram_plot_data]=...
    variogram_estimation(z_tilde,lon,lat,binendpoints,V_eps);
% save variogram_estimates sig2_eps_variogram ...
%     sig2_xi_variogram variogram_plot_data
% load variogram_estimates


%%%% other MM estimation

% load variogram_estimates sig2_xi_variogram

data_MM=[lat lon z_tilde];
K_MM=binest(data_MM,So,sig2_eps,V_eps,sig2_xi_variogram);

% save MMestimates K_MM sig2_xi_variogram


%%%% EM estimation

[K_em sig2_xi_em T]=EM(So,z_tilde,V_eps,sig2_eps);

% save emestimates K_em sig2_xi_em


%%%%%%%%%   obtain predictions  %%%%%%%%%

data=[lat lon z_tilde];
pred_locs=[pred_lat pred_lon];

% MM
% load MMestimates
[pred_notrend_MM sig2_MM]=FRK(data,pred_locs,K_MM,sig2_xi_variogram,...
    sig2_eps,V_eps,1);

% EM
% load emestimates
[pred_notrend_EM sig2_EM]=FRK(data,pred_locs,K_em,sig2_xi_em,...
    sig2_eps,V_eps,1);

% add back trend
Xp=[ones(length(pred_lon),1) pred_lon];
pred_MM=pred_notrend_MM+Xp*betahat;
pred_EM=pred_notrend_EM+Xp*betahat;

% save predictions pred_MM sig2_MM pred_EM sig2_EM



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   plots for the document  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bord_proj=dlmread('proj_borders.csv')*90;
longs=[-180:180 ones(1,179)*179.99 180:(-1):-180 -ones(1,179)*180];
lats=[-ones(1,359)*89.99 -89.99 -89:89 89.99 ...
    ones(1,359)*89.99 89.99 89:(-1):-89 -89.99];
[lonhull lathull]=mollweide(longs',lats');
gi=.3;    
temp=jet(101); offse=7;
noeljet=temp((offse+1):(end-offse),:);


%%% the data

markersizes=90-abs(lat);
[x y]=mollweide(lon,lat);
scatter(x, y,markersizes,z,'.')
    xlim([-181 181]); ylim([-91 91])
    title('Data','FontSize',14)
hold on, 
    plot(bord_proj(:,1),bord_proj(:,2),'-','Color',[gi gi gi])
    plot(lonhull,lathull,'-','Color',[gi gi gi])        
hold off
axis off
colorbar; colormap(noeljet)
set(gcf, 'color', 'white')
% delete('../plots/data.png')
% export_fig('../plots/data.png','-painters','-m1')


%%% examine normality

subplot(1,3,1), hist(z_tilde), title('Histogram of Detrended Data')
subplot(1,3,2), kde(z_tilde), 
    title('Smoothed Histogram of Detrended Data')
subplot(1,3,3), qqplot(z_tilde), title('Normal Q-Q Plot of Detrended Data')    
    ylabel('Quantiles of Detrended Data')
set(gcf, 'color', 'white')
% delete '../plots/normality.pdf'; export_fig '../plots/normality' '-pdf'


%%% bisquare basis functions

% one two-dimensional bisquare function
length=49;
mid=(length+1)/2;
radius=mid;
bisquare2d=zeros(length,length);
ind=1;
for i=1:length
    for j=1:length
        if norm([i; j]-[mid; mid])<=radius
            bisquare2d(i,j)=(1-(norm([i; j]-[mid; mid])/radius)^2)^2 ;
        end
    end
end
subplot(1,2,1)
    surf(1:length,1:length,bisquare2d)
    axis tight, set(gca,'XTickLabel','','YTickLabel','')
subplot(1,2,2)
    imagesc(bisquare2d), colorbar
    set(gca,'XTickLabel','','YTickLabel','')
set(gcf, 'color', 'white')
% delete '../plots/bisquare2d.pdf'; export_fig '../plots/bisquare2d' '-pdf'
clear length

% several one-dimensional bisquare functions
center={128.5,[32.5 96.5 160.5 224.5]};
radius=1.5*[96 32]; locs=1:256;
S=zeros(256,5);
i=1;
for l=1:length(center),
    for jl=1:length(center{l}),
        for u=1:256,
            if abs(locs(u)-center{l}(jl))<=radius(l)
                S(u,i)=(1-(abs(locs(u)-center{l}(jl))/radius(l))^2)^2 ;
            end
        end
        i=i+1;
    end
end
subplot(1,2,1), plot(S), axis tight
subplot(1,2,2), plot(S*normrnd(0,1,5,1)), axis tight
set(gcf, 'color', 'white')
% delete '../plots/bisquare1d.pdf'; export_fig '../plots/bisquare1d' '-pdf'


%%% the basis function centers for the XCO2 example

gi=.75;
wi=2;
bord_proj=dlmread('proj_borders.csv')*90;
longs=[-180:180 ones(1,179)*179.99 180:(-1):-180 -ones(1,179)*180];
lats=[-ones(1,359)*89.99 -89.99 -89:89 89.99 ...
    ones(1,359)*89.99 89.99 89:(-1):-89 -89.99];
[lonhull lathull]=mollweide(longs',lats');

[x y]=mollweide(level1(:,1),level1(:,2));
scatter(x, y, 230,'r*','LineWidth',wi);  
hold on,
    [x y]=mollweide(level2(:,1),level2(:,2));
    scatter(x, y, 130,'bx','LineWidth',wi);
    [x y]=mollweide(level3(:,1),level3(:,2));
    scatter(x, y, 50,'g+','LineWidth',wi);    
hold off
legend('Res 1','Res 2','Res 3')

xlim([-182 181]); ylim([-92 91]);    
hold on, 
    plot(bord_proj(:,1),bord_proj(:,2),'-','Color',[gi gi gi])
    plot(lonhull,lathull,'-k')
    [x y]=mollweide(level3(:,1),level3(:,2));
    scatter(x, y, 50,'g+','LineWidth',wi);    
    [x y]=mollweide(level2(:,1),level2(:,2));
    scatter(x, y, 130,'bx','LineWidth',wi);
    [x y]=mollweide(level1(:,1),level1(:,2));
    scatter(x, y, 230,'r*','LineWidth',wi);
hold off
axis off
set(gcf, 'color', 'white')
% delete 'plots/basisfunctioncenters.pdf'
% export_fig 'plots/basisfunctioncenters' '-pdf'


%%% the semivariogram

% variogram_plot_data={centers,rob_meanb,beta_variogram,n_pairsb,bins};
centers=variogram_plot_data{1};
rob_semivar=variogram_plot_data{2};
beta_variogram=variogram_plot_data{3};
n_pairs=variogram_plot_data{4};
bins=variogram_plot_data{5};
plot(centers,rob_semivar,'xr')
    xlim([0 bins(end)]), ylim([0 max(rob_semivar)*1.05])
    hold on
        plot([0 centers],beta_variogram(1)+beta_variogram(2)*[0 centers])
    hold off
    for i=1:length(n_pairs),
        text(centers(i)-1.5,rob_semivar(i)*.94,...
            strcat('(',num2str(round(n_pairs(i))),')'),...
            'FontSize',14)
    end
    xlabel('great-arc distance (in km)')
set(gcf, 'color', 'white')    
% delete '../plots/semivariogram2.pdf'
% export_fig '../plots/semivariogram2' '-pdf'


%%% the estimates of the matrix K

temp=flipud(lbmap(101,'BlueRed'));
noelbluered=[temp(10:1:40,:);temp(45,:);ones(1,3);temp(57,:);...
    temp(62:1:92,:)];

ma_mm=max(abs(vec(K_MM)));
ma_em=max(abs(vec(K_em)));

subplot(1,1,1), imagesc(K_MM), colorbar
    title('MM estimate of K','FontSize',14)
    vline(32.5,'k-'), hline(32.5,'k-')
    vline(124.5,'k-'), hline(124.5,'k-') 
set(gcf, 'color', 'white')
caxis([-ma_mm ma_mm]), colormap(noelbluered)
% delete '../plots/KMM.pdf'; export_fig '../plots/KMM.pdf'

subplot(1,1,1), imagesc(K_em), colorbar
    title('EM estimate of K','FontSize',14)
    vline(32.5,'k-'), hline(32.5,'k-')
    vline(124.5,'k-'), hline(124.5,'k-') 
set(gcf, 'color', 'white')
caxis([-ma_em ma_em]), colormap(noelbluered)
% delete '../plots/KEM.pdf'; export_fig '../plots/KEM.pdf'


% comparison of the elements
r=size(K_em,1);
bet=regress(vec(K_MM),[ones(numel(K_em),1) vec(K_em)]);
scatter(vec(K_em),vec(K_MM),25,vec(eye(r)),'.')
axis tight, xlabel('K_{EM}'), ylabel('K_{MM}')
hold on,
    scatter(diag(K_em),diag(K_MM),25,'r.')
    plot([min(K_MM) max(K_MM)],[min(K_MM) max(K_MM)],'m')
    plot([min(K_MM) max(K_MM)],bet(1)+[min(K_MM) max(K_MM)]*bet(2),'g')
    vline(0,'k-'), hline(0,'k-')
hold off
legend({'off-diag','diag','x=y','OLS line'},'Location','NorthWest')
set(gcf, 'color', 'white')
% delete('../plots/Kelements.png')
% export_fig('../plots/Kelements.png','-painters','-m2')

% comparison of the estimated correlation matrices
Rmm=diag(diag(K_MM))^(-.5)*K_MM*diag(diag(K_MM))^(-.5);
Rem=diag(diag(K_em))^(-.5)*K_em*diag(diag(K_em))^(-.5);
r=size(Rem,1);
bet=regress(vec(Rmm),[ones(numel(Rem),1) vec(Rem)]);
scatter(vec(Rem),vec(Rmm),25,vec(eye(r)),'.')
axis tight, xlabel('R_{EM}'), ylabel('R_{MM}')
hold on,
    scatter(diag(Rem),diag(Rmm),25,'r.')
    plot([min(Rmm) max(Rmm)],[min(Rmm) max(Rmm)],'m')
    plot([min(Rmm) max(Rmm)],bet(1)+[min(Rmm) max(Rmm)]*bet(2),'g')
    vline(0,'k-'), hline(0,'k-')
hold off
legend({'off-diag','diag','x=y','OLS line'},'Location','NorthWest')
set(gcf, 'color', 'white')
% delete('../plots/correlements.png')
% export_fig('../plots/correlements.png','-painters','-m2')


%%% the predictions and standard errors

% load predictions
[pred_x pred_y]=mollweide(pred_lon,pred_lat);
markersizes_pred=90-abs(pred_lat);
[x y]=mollweide(lon,lat);

mipred=min([z; pred_MM]); mapred=max([z; pred_MM]);
mise=min(sig2_MM)^.5; mase=max(sig2_MM)^.5;

% for MM 

ax(1)=subplot(3,1,1);
scatter(x, y,markersizes,z,'.')
    xlim([-181 181]); ylim([-91 91])
    title('Data','FontSize',14)
hold on, 
    plot(bord_proj(:,1),bord_proj(:,2),'-','Color',[gi gi gi])
    plot(lonhull,lathull,'-','Color',[gi gi gi])        
hold off
axis off

ax(2)=subplot(3,1,2);
scatter(pred_x, pred_y,markersizes_pred,pred_MM,'.')
    xlim([-181 181]); ylim([-91 91])
    title('FRK predictions using MM estimates','FontSize',14)
hold on, 
    plot(bord_proj(:,1),bord_proj(:,2),'-','Color',[gi gi gi])
    plot(lonhull,lathull,'-','Color',[gi gi gi])        
hold off
axis off;  colorbar
h1=colorbar; caxis([mipred mapred])

ax(3)=subplot(3,1,3);
scatter(pred_x, pred_y,markersizes_pred,sig2_MM.^.5'.')
    xlim([-181 181]); ylim([-91 91]), caxis([mise mase])
    title('FRK standard errors using MM estimates','FontSize',14)
hold on, 
    plot(bord_proj(:,1),bord_proj(:,2),'-','Color',[gi gi gi])
    plot(lonhull,lathull,'-','Color',[gi gi gi])        
hold off
axis off;  colorbar
set(gcf, 'color', 'white')
h2=colorbar;
        
set(h1, 'Position', [.925 0.42 .03 .52])
set(h2, 'Position', [.925 0.12 .03 .22])
colormap(noeljet)
set(gcf, 'color', 'white')
% delete('../plots/MMpred2.png')
% export_fig('../plots/MMpred2.png','-painters','-m1.5')


% for EM 

mipred=min([z; pred_EM]); mapred=max([z; pred_EM]);
mise=min(sig2_EM)^.5; mase=max(sig2_EM)^.5;

ax(1)=subplot(3,1,1);
scatter(x, y,markersizes,z,'.')
    xlim([-181 181]); ylim([-91 91])
    title('Data','FontSize',14)
hold on, 
    plot(bord_proj(:,1),bord_proj(:,2),'-','Color',[gi gi gi])
    plot(lonhull,lathull,'-','Color',[gi gi gi])        
hold off
axis off

ax(2)=subplot(3,1,2);
scatter(pred_x, pred_y,markersizes_pred,pred_EM,'.')
    xlim([-181 181]); ylim([-91 91])
    title('FRK predictions using EM estimates','FontSize',14)
hold on, 
    plot(bord_proj(:,1),bord_proj(:,2),'-','Color',[gi gi gi])
    plot(lonhull,lathull,'-','Color',[gi gi gi])        
hold off
axis off;  colorbar
h1=colorbar; caxis([mipred mapred])

ax(3)=subplot(3,1,3);
scatter(pred_x, pred_y,markersizes_pred,sig2_EM.^.5'.')
    xlim([-181 181]); ylim([-91 91]), caxis([mise mase])
    title('FRK standard errors using EM estimates','FontSize',14)
hold on, 
    plot(bord_proj(:,1),bord_proj(:,2),'-','Color',[gi gi gi])
    plot(lonhull,lathull,'-','Color',[gi gi gi])        
hold off
axis off;  colorbar
set(gcf, 'color', 'white')
h2=colorbar;
        
set(h1, 'Position', [.925 0.42 .03 .52])
set(h2, 'Position', [.925 0.12 .03 .22])
colormap(noeljet)
set(gcf, 'color', 'white')
% delete('../plots/EMpred.png')
% export_fig('../plots/EMpred.png','-painters','-m1.5')




%%% compare empirical to (estimated) theoretical semivariograms

% load MMestimates K_MM sig2_xi_variogram
% load emestimates K_em sig2_xi_em

level{1}=level1; level{2}=level2; level{3}=level3;
for i=1:3, [B,IX] = sort(level{i},1); level{i}=level{i}(IX(:,2),:); end
cents=[level{1}; level{2}; level{3}];
radii=[repmat(6241.10994046,32,1); repmat(3491.03461571,92,1); ...
    repmat(2047.58161958,272,1)];

Sp=bisquare(distance_spherical([pred_lon pred_lat],cents),radii');

refpoint=[100, 30];
lonbins=((-5):10:55) + refpoint(1);
latbins=((-5):10:55) + refpoint(2);

ns=zeros(length(lonbins)-1,4);
semiv_lon_emp=zeros(length(lonbins)-1,1);
semiv_lat_emp=zeros(length(lonbins)-1,1);
semiv_lon_EM=zeros(length(lonbins)-1,1);
semiv_lat_EM=zeros(length(lonbins)-1,1);
semiv_lon_MM=zeros(length(lonbins)-1,1);
semiv_lat_MM=zeros(length(lonbins)-1,1);

loncen=mean([lonbins(2:end);lonbins(1:(end-1))]);
latcen=mean([latbins(2:end);latbins(1:(end-1))]);

refdata=z(lon>lonbins(1)&lon<lonbins(2) & lat>latbins(1)&lat<latbins(2));
b_ref=bisquare(distance_spherical([loncen(1) latcen(1)],cents),radii');

for i=1:(length(lonbins)-1),
    
    % lon
    compdata=z(lon>lonbins(i) & lon<lonbins(i+1) & ...
        lat>latbins(1) & lat<latbins(2));
    diff=vec(pdist2(refdata,compdata)); ns(i,1)=length(diff);
    semiv_lon_emp(i)=(.5*mean(sqrt(diff))^4/(.457+.494/length(diff))).^.5;
    
    b_comp=bisquare(distance_spherical([loncen(i) latcen(1)],cents),...
        radii'); b_diff=b_comp-b_ref;
    semiv_lon_EM(i)=(.5*(b_diff*K_em*b_diff'+2*sig2_xi_em+2*sig2_eps)).^.5;

    semiv_lon_MM(i)=(.5*(b_diff*K_MM*b_diff'+...
        2*sig2_xi_variogram+2*sig2_eps)).^.5;
    
    % lat
    compdata=z(lon>lonbins(1) & lon<lonbins(2) & ...
        lat>latbins(i) & lat<latbins(i+1));
    diff=vec(pdist2(refdata,compdata)); ns(i,2)=length(diff);
    semiv_lat_emp(i)=(.5*mean(sqrt(diff))^4/(.457+.494/length(diff))).^.5;
    
    b_comp=bisquare(distance_spherical([loncen(1) latcen(i)],cents),...
        radii'); b_diff=b_comp-b_ref;
    semiv_lat_EM(i)=(.5*( (latcen(i)-latcen(1))*betahat(2)+...    
        b_diff*K_em*b_diff'+2*sig2_xi_em+2*sig2_eps)).^.5;
    semiv_lat_MM(i)=(.5*( (latcen(i)-latcen(1))*betahat(2)+...    
        b_diff*K_MM*b_diff'+2*sig2_xi_variogram+2*sig2_eps)).^.5;
    
end

% plot
ma=max([semiv_lon_emp(:); semiv_lat_emp(:); semiv_lon_EM(:);...
    semiv_lat_EM(:); semiv_lon_MM(:); semiv_lat_MM(:)]);
subplot(1,2,1), plot(loncen,semiv_lon_emp,'o'), xlabel('lon')
    % ylim([0 ma])
    hold on, 
        plot(loncen,semiv_lon_MM,'r--')
        plot(loncen,semiv_lon_EM,'k-')
    hold off
    xlabel('longitude')
%         legend({'Data','MM','EM'},'Location','SouthEast')
subplot(1,2,2), plot(latcen,semiv_lat_emp,'o'), xlabel('lat')
    % ylim([0 ma])
    hold on, 
        plot(latcen,semiv_lat_EM,'k-')
        plot(latcen,semiv_lat_MM,'r--')
    hold off
    xlabel('latitude')
%         legend({'Data','MM','EM'},'Location','SouthEast')
   
set(gcf, 'color', 'white')
% delete(['../plots/corrcomp_ref_',num2str(refpoint),'.pdf'])
% export_fig(['../plots/corrcomp_ref_',num2str(refpoint),'.pdf'])

