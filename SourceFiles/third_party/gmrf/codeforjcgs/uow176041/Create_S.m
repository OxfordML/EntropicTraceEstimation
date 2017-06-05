function S=Create_S(crd,level1,level2,level3,type)

if nargin<4, type=1; end % type=1 is bisquare, type=2 is Poisson

level{1}=level1;
level{2}=level2;
level{3}=level3;

for i=1:3
    [B,IX] = sort(level{i},1);
    level{i}=level{i}(IX(:,2),:);
end


% MAKING THE S MATRIX

if type==1, % bisquare
    
    % radii agreed upon with hai:
    rl=[6241.10994046; 3491.03461571; 2047.58161958];

    S=zeros(size(crd,1),size(level1,1)+size(level2,1)+size(level3,1));
    count=0;
    for i=1:3
        % hrl=distance_spherical(level{i},level{i});
        % rl(i,1)=1.5*min(hrl(hrl>1e-3));
        for j=1:length(level{i})
            count=count+1;
            h=distance_spherical(level{i}(j,:),crd);
            s=(1-(h./rl(i,1)).^2).^2;
            s(h>rl(i,1))=0;
            S(:,count)=s;       
        end
    end
    
end

if type==2, % Poisson kernel
    
    earth_radius=6371.0087714; % WGS84 mean radius
    eff_range=.25;

    S=zeros(size(crd,1),size(level1,1)+size(level2,1)+size(level3,1));
    count=0;
    for i=1:3
        hrl=distance_spherical(level{i},level{i});
        r_bisquare=1.5*min(hrl(hrl>1e-3));
        a=1-cos(r_bisquare/earth_radius)*eff_range^(2/3);
        b=1-eff_range^(2/3);    
        r(i,1)=(a-sqrt(a^2-b^2))/b;
        for j=1:length(level{i})
            count=count+1;
            x=cos(distance_spherical(level{i}(j,:),crd)/earth_radius);
            s=(1-r(i,1))^2*(1+r(i,1))^(-1)*(1-r(i,1)^2)./(1-2.*x.*r(i,1)+r(i,1)^2).^1.5;
            S(:,count)=s;       
        end
    end
    
end
