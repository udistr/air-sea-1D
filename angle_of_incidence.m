function y=angle_of_incidence(lat,lon,jd,time)

delta=(asin(-sin(23.45./180.*pi).*cos(360./365.25.*(jd+10)./180.*pi)));
lamda=lat./180.*pi;
omega=time./24.*2.*pi+lon./360*2*pi+pi;

y=cos(delta).*cos(lamda).*cos(omega)+sin(delta).*sin(lamda);