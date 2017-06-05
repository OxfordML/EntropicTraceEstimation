% function to bin observations on the globe

function [W]=bin(lon,lat)

W={};
    
%%%%%%%%%   6x5 cells per bin   %%%%%%%%%%%%%%%%%%%%    
    
cent_lon=[]; cent_lat=[]; count=[];
i=1;

for centlons=-172.5:7.5:172.5,
    
    % regular bins
    for centlats=-75:5:75,
        cent_lon(i,1)=centlons;
        cent_lat(i,1)=centlats;
        W{i}=find(abs(lon-centlons)<3.75 & abs(lat-centlats)<2.5 );
        count(i,1)=length(W{i});
        i=i+1;
    end
    
    % pole bins
    cent_lon(i,1)=centlons; cent_lat(i,1)=83.75; 
        W{i}=find(abs(lon-centlons)<3.75 & lat>77.5); count(i,1)=length(W{i}); i=i+1;
    cent_lon(i,1)=centlons; cent_lat(i,1)=-83.75; 
        W{i}=find(abs(lon-centlons)<3.75 & lat<-77.5); count(i,1)=length(W{i}); i=i+1;
        
end

% bins centered opposite greenwich
for centlats=-75:5:75,
    cent_lon(i,1)=180;
    cent_lat(i,1)=centlats;
    W{i}=find(abs(lon)>176.25 & abs(lat-centlats)<2.5);
    count(i,1)=length(W{i});
    i=i+1;
end
    cent_lon(i,1)=180; cent_lat(i,1)=83.75; 
        W{i}=find(abs(lon)>176.25 & lat>77.5); count(i,1)=length(W{i}); i=i+1;
    cent_lon(i,1)=180; cent_lat(i,1)=-83.75; 
        W{i}=find(abs(lon)>176.25 & lat<-77.5); count(i,1)=length(W{i});    
    
        
% discard all bins with less than 15 observations        
W=W(count>=15);



%%%%%%%%%   combine bins to achieve >=15 obs per bin   %%%%%%%%%%%%%%%%%%%%    


% while (sum(count<15)>0),
%     
%     lons=unique(cent_lon);
%     lats=unique(cent_lat);
% 
%     for i=1:length(lats),
%         for j=1:length(lons),
%         
%             index=find(cent_lon==lons(j) & cent_lat==lats(i));
%         
%             if count(index)<15,
%             
%                 % index of nearest neighbor
%                 dists=distance_spherical([lons(j) lats(i)],[cent_lon cent_lat]);
%                 [sortdist sortind]=sort(dists);
%                 nextind=sortind(2);
%                 nextindlon=cent_lon(nextind); nextindlat=cent_lat(nextind);
% 
%                 % delete entries of bin and neighbor bin
%                 cent_lon([index nextind])=[];
%                 cent_lat([index nextind])=[];
%                 w_temp=[W{index}; W{nextind}];
%                 W([index nextind])=[];
%                 count([index nextind])=[];
% 
%                 % add combined bin at the bottom
%                 cent_lon=[cent_lon; mean([lons(j) nextindlon])];
%                 cent_lat=[cent_lat; mean([lats(i) nextindlat])];
%                 W=[W {w_temp}];
%                 count=[count; length(w_temp)];          
% 
%             end
%         end
%     end
% end
