function cord = LatLong2MCMF(lat, lon)
%{
This function gets MCMF cordinate from moon latitude and longitude.
Lat = latitude in degrees
Lon = longitude in degrees

cord = [x,y,z] cordinates in Moon centered Moon fixed frame.
%}
Rm = 1740; %Radius of moon in km
z = Rm*sind(lat);
x = Rm*cosd(lat).*cosd(lon);
y = Rm*cosd(lat).*sind(lon);
cord = [x,y,z];%R vector in ECEF cord in kms
end