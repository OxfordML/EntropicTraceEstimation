%%% Mollweide projection of long/lat coordinates
% author: matthias katzfuss

function [x y]=mollweide(lon,lat)

lambda=lon*pi/180;
phi=lat*pi/180;


% find auxiliary angle theta
eps=1e-5;
thetap_old=2*asin(2*phi/pi);
done=0;
    
while done==0,
    thetap_new=thetap_old-(thetap_old+sin(thetap_old)-pi*sin(phi))./...
        (1+cos(thetap_old));
    if max(abs(thetap_new-thetap_old))<eps, done=1; end
    thetap_old=thetap_new;
end

theta=.5*thetap_new;


% calculate the projected coordinates
x=2*sqrt(2)*lambda.*cos(theta)*(90/4.427456261289762);
y=sqrt(2)*sin(theta)*(90/1.413651161814046);
