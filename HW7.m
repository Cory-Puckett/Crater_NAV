%{
Homework 7/Project main file

Assumptions
- pixels are square
- no shearing
- no distortion
- camera pointing is center of image

estimated altitude of 363 km based on napkin math
%}
clearvars
close all
load moon_Image.mat
Rm = 1740; %Radius of moon in km

%image size
imgWidth = 4034;
imgHeight = 3026;

up = imgWidth/2;
vp = imgHeight/2;
uFOV = deg2rad(20); %Horizontal FOV given by assignment
dx = imgWidth/(2*tan(uFOV/2));
dy = dx;
K = [dx, 0, up;
    0, dy, vp;
    0, 0, 1]; % Camera calibration matrix

% Crater database reduction, but don't actually need
% maxCraterDia = 55; %Maximum crater diameter to keep and look for
% minCraterDia = 2; %Minimum crater diameter to keep and look for
% craterDatabase =  DataReduce(maxCraterDia,minCraterDia,5,10);

%Image cordinates in upside down image
upsidedownCord = [1353,1849;191,852;3417,1534;1561,1494];
imgCord = [(imgWidth-upsidedownCord(:,1)),(imgHeight-upsidedownCord(:,2)), ...
    ones(length(upsidedownCord),1)];

craterLatLon = [63.633,93.6;62.966,89.980;65.694,95.7;63.975,93.087];
% craterLatLon = [63.633,93.6;62.966,89.980;65.694,95.7;63.975,93.0];
%need lat long to mcmf func

% Plotting of points used.
imshow(moon_img);
hold on;
scatter(imgCord(:,1),imgCord(:,2),"filled",DisplayName="Crater points used");
legend;
camroll(-180); %Flipping image so that craters look like craters instead of mountains.
exportgraphics(gca,"Crater Points.png",Resolution=600);

%% Getting line of sight vectors
%LOS vectors with cam x, y, z as columns
u = zeros(size(imgCord)); %to be LOS unit vectors
%Looping through each point to compute LOS.
for i = 1:size(imgCord,1)
    u_i = K^-1*imgCord(i,:)';
    u(i,:) = u_i'/norm(u_i);
end

%% Getting MCMF cords of craters
%World points picked x, y, z as columns
P = LatLong2MCMF(craterLatLon(:,1),craterLatLon(:,2));

%% Grunert's method
[TMCMF2cam, r] = Grunerts(u,P);
r2LatLon(r, Rm);
% r = [2.027569592758852e+02;7.540225095545408e+02;1.759188745910317e+03];
% TMCMF2cam = [0.319751603053375,0.848624126748070,-0.421421408860385;0.529644580363958,-0.528878548310257,-0.663147117635538;-0.785643386553902,-0.011161211328058,-0.618578933140958];
u_exp = zeros(size(imgCord)); %to be expected LOS unit vectors
imgCord_exp = zeros(size(imgCord));
for i = 1:length(u)
    ui_exp = ((TMCMF2cam*(P(i,:)'-r))/norm(TMCMF2cam*(P(i,:)'-r)));
    u_exp(i,:) = ui_exp';
    imgCord_exp(i,:) = (K*(ui_exp/ui_exp(3)))';
end

%% Plotting with expected image cordinates
figure;
imshow(moon_img);
hold on;
scatter(imgCord(:,1),imgCord(:,2),"filled",DisplayName="Crater points used");
scatter(imgCord_exp(:,1),imgCord_exp(:,2),"filled",DisplayName="Expected crater points");
legend;
camroll(-180); %Flipping image so that craters look like craters instead of mountains.
exportgraphics(gca,"Result.png",Resolution=600);


%% Plotting craters in area
load craterDataBaseSmall.mat craterDatabase

craterDatabase;






