function d=delta(n)

%Declination Angle
%This function gets the Julian day and calculate declination of the sun

d=(asin(-sin(23.45./180.*pi).*cos(360./365.25.*(n+10)./180.*pi)));


end
