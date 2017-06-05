load('uow176043.mat');
pred_lat(pred_lat == 89.75) = 90;
pred_lat(pred_lat == -89.75) = -90; 
lat2 = (lat - min(pred_lat))/1 + 1;
lon2 = (lon - min(pred_lon))/1.25 + 1;
obsMap2 = zeros(288,181);
for i = 1:26633
obsMap2(lon2(i),lat2(i)) = 1;
end;

subdata2 = zeros(288,181);
for i = 1:26633
subdata2(lon2(i),lat2(i)) = z(i);
end;