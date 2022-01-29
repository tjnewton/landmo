
%% coastsurf.m
% generates interpolated vertical land velocity and estimates uncertainty
% required files:
% washington_coastline.mat           - PNW coastline file
% POI.xlsx               - points of interest for histograms
% profiledata.xlsx       - list of points for velocity profiles
% new.xlsx       - complete VLM dataset
% WAcoastline2dec.xlsx   - coastline points for interpolation
% WAmodelpointsRay6.xlsx - model points

% gridded model point generation
close all
clear all

[num] = xlsread('WAmodelpointsRay6.xlsx'); % reads data file and imports all elements that are numbers
latitude = num(:, 2); 
longitude = num(:, 1);
velocity = num(:, 3);

F = scatteredInterpolant(longitude, latitude, velocity, 'linear');
[xq,yq] = meshgrid(-125:0.1:-121.5, 45.5:0.1:49.0); % creates query points for interpolated plot
vq = F(xq, yq);
figure
surf(xq, yq, vq, 'EdgeColor','none')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('Tectonic Model')
colormap jet
xlim([-125 -121.5])
ylim([45.5 49])
zlim([-2 4])
colorbar
caxis([-3 3])

modelx =[-124.95:0.3:-121.35];
modelx = modelx.';
modely = [45.45:0.2:49.05];
modely = modely.';
catarray = [];

for i = 1:13
    data = [];
    data(1:19,1) = modelx(i);
    data(1:19,2) = modely;
    catarray = vertcat(catarray, data);
end

interpvelarray = [];
for i=1:length(catarray(:,1))
[rowx,colx] = find(abs(abs(xq-catarray(i,1)) == min(abs(xq(:)-catarray(i,1)))));
[rowy,coly] = find(abs(abs(yq-catarray(i,2)) == min(abs(yq(:)-catarray(i,2)))));
interpVel = vq(rowy(1),colx(1));
interpvelarray = vertcat(interpvelarray, interpVel);
end
catarray(:,3) = interpvelarray;
% xlswrite('modeldata.xlsx',catarray);


%% vertical velocity interpolation 
close all
clear all

[num, txt, raw] = xlsread('new.xlsx'); % reads data file and imports all elements that are numbers and text
numNaN = any(isnan(num), 1); % finds NaN values in 'num'
num(:, numNaN) = []; % deletes columns with NaN values in 'num'
latlat = num(:, 2); 
lonlon = num(:, 1);
velvel = num(:, 3);
typtyp = txt(:, 6);

% exclude data of quality index 2 and 4
exclude2 = find(num(:,5) == 2);
num(exclude2,:) = [];
txt(exclude2,:) = [];
exclude4 = find(num(:,5) == 4);
num(exclude4,:) = [];
txt(exclude4,:) = [];

latitude = num(:, 2); 
longitude = num(:, 1);
velocity = num(:, 3);
standardDeviation = num(:, 4);
qualityindex = num(:,5);
txt(1,:) = [];
station = txt(:, 1); 
type = txt(:, 6);

typeref = type;
longituderef = longitude;
latituderef = latitude;
velocityref = velocity;
uncertref = standardDeviation;

[num, txt] = xlsread('modeldata.xlsx'); % reads model data
modellatitude = num(:, 2); 
modellongitude = num(:, 1);
modelvelocity = num(:, 3);
modeltype = txt(:,1);
latitude = vertcat(latitude, modellatitude); 
longitude = vertcat(longitude, modellongitude); 
velocity = vertcat(velocity, modelvelocity); 
type = vertcat(type, modeltype); 

F = scatteredInterpolant(longitude, latitude, velocity); % creates interpolation function
[xq,yq] = meshgrid(-126:0.1:-121, 45:0.1:49.5); % creates query points for interpolated plot
vq = F(xq, yq);

% typeref = type;
% longituderef = longitude;
% latituderef = latitude;
% velocityref = velocity;
% uncertref = standardDeviation;
%end filtering

% removing outliers for GIA gradient figure
checkme(:,1) = longitude;
checkme(:,2) = velocity;
findy = find(checkme(:,1) > -123 & checkme(:,2) < -1);
longitude(findy) = [];
latitude(findy) = [];
velocity(findy) = [];

F = scatteredInterpolant(longitude, latitude, velocity);
[xq,yq] = meshgrid(-126:0.1:-121, 45:0.1:49.5);
vq = F(xq, yq);
figure
surf(xq, yq, vq, 'EdgeColor','none') % interpolated vertical velocity surface
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
%title('Interpolated Velocity')
title('Interpolated Vertical Velocity')
xlim([-125 -121.5])
ylim([45.5 49])
zlim([-3 4])

contour20 = load('cas_slab2_dep_20.txt');
contour20(:,1) = (contour20(:,1) - 360);
zz20 = zeros(1,length(contour20));
zz20 = zz20+4;
%plot3(contour20(:,1),contour20(:,2), zz20,'k', 'LineWidth', 2, 'DisplayName','20 km contour','LineStyle','-')

contour40 = load('cas_slab2_dep_40.txt');
contour40(:,1) = (contour40(:,1) - 360);
zz40 = zeros(1,length(contour40));
zz40 = zz40+4;
%plot3(contour40(:,1),contour40(:,2), zz40,'k', 'LineWidth', 2, 'DisplayName','40 km contour','LineStyle','--')

contour60 = load('cas_slab2_dep_60.txt');
contour60(:,1) = (contour60(:,1) - 360);
zz60 = zeros(1,length(contour60));
zz60 = zz60+4;
%plot3(contour60(:,1),contour60(:,2), zz60,'k', 'LineWidth', 2, 'DisplayName','60 km contour','LineStyle',':')

contour80 = load('cas_slab2_dep_80.txt');
contour80(:,1) = (contour80(:,1) - 360);
zz80 = zeros(1,length(contour80));
zz80 = zz80+4;
%plot3(contour80(:,1),contour80(:,2), zz80,'k', 'LineWidth', 2, 'DisplayName','80 km contour','LineStyle','-.')

contour100 = load('cas_slab2_dep_100.txt');
contour100(:,1) = (contour100(:,1) - 360);
zz100 = zeros(1,length(contour100));
zz100 = zz100+4;
%plot3(contour100(:,1),contour100(:,2), zz100,'k', 'LineWidth', 2, 'DisplayName','100 km contour','LineStyle','-')

load PNWcoast.dat;
zz = zeros(1,length(PNWcoast));
zz = zz+4;
plot3(PNWcoast(:,1),PNWcoast(:,2), zz,'k', 'LineWidth', 1, 'DisplayName','Coastline')
h = colorbar;
caxis([-3 3])
ylabel(h, 'Vertical Velocity (mm/yr)')
pbaspect([1 1 0.3])

% to export .asc file
%arcgridwrite('arc_vel_surf_new.asc',xq(1,:),yq(:,1),vq);
%arcgridwrite('interpolated_vertical_velocity_0.01.asc',xq(1,:),yq(:,1),vq);

%%

% plot with all data

[num, txt] = xlsread('new.xlsx'); % reads data file and imports all elements that are numbers
numNaN = any(isnan(num), 1); % finds NaN values in 'num'
num(:, numNaN) = []; % deletes columns with NaN values in 'num'

% exclude data of quality index 2 and 4
exclude2 = find(num(:,5) == 2);
num(exclude2,:) = [];
txt(exclude2,:) = [];
exclude4 = find(num(:,5) == 4);
num(exclude4,:) = [];
txt(exclude4,:) = [];

latitude2 = num(:, 2); 
longitude2 = num(:, 1);
velocity2 = num(:, 3);
latitude6 = num(:, 2); 
longitude6 = num(:, 1);
velocity6 = num(:, 3);
type6 = txt(:, 6);

latlonvel(:,1) = latitude6;
latlonvel(:,2) = longitude6;
latlonvel(:,3) = velocity6;


[num] = xlsread('modeldata.xlsx');
modellatitude = num(:, 2); 
modellongitude = num(:, 1);
modelvelocity = num(:, 3);
modellatitude = modellatitude;
modellongitude = modellongitude;
latitude2 = vertcat(latitude2, modellatitude); 
longitude2 = vertcat(longitude2, modellongitude); 
velocity2 = vertcat(velocity2, modelvelocity);

index = find(contains(type,'GPS'));
GPSlatitude = latitude2(index);
GPSlongitude = longitude2(index);
GPSvelocity = velocity2(index);
index = find(contains(type,'Leveling'));
LEVlatitude = latitude2(index);
LEVlongitude = longitude2(index);
LEVvelocity = velocity2(index);
index = find(contains(type,'WL'));
WLlatitude = latitude2(index);
WLlongitude = longitude2(index);
WLvelocity = velocity2(index);
index = find(contains(type,'MOD'));
MODlatitude = latitude2(index);
MODlongitude = longitude2(index);
MODvelocity = velocity2(index);

z = zeros(length(latitude2),1);
z = z+5;
z1 = zeros(length(GPSlatitude),1);
z1 = z1+5;
z2 = zeros(length(LEVlatitude),1);
z2 = z2+5;
z3 = zeros(length(WLlatitude),1);
z3 = z3+5;
z4 = zeros(length(MODlatitude),1);
z4 = z4+5;

F = scatteredInterpolant(longitude, latitude, velocity);
[xq,yq] = meshgrid(-126:0.1:-121, 45:0.1:49.5); % creates query points for interpolated plot
vq = F(xq, yq);
figure
surf(xq, yq, vq, 'EdgeColor','none','DisplayName','Interpolated Surface')
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('Interpolated Velocity')
scatter3(MODlongitude, MODlatitude, z4, 40, MODvelocity,'o', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Model','LineWidth',0.1,'MarkerEdgeAlpha',0.6)
scatter3(GPSlongitude, GPSlatitude, z1, 70, GPSvelocity,'<', 'filled', 'MarkerEdgeColor', 'k','DisplayName','GPS','LineWidth',0.1,'MarkerEdgeAlpha',0.6)
scatter3(LEVlongitude, LEVlatitude, z2, 70, LEVvelocity,'s', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Leveling','LineWidth',0.1,'MarkerEdgeAlpha',0.6)
scatter3(WLlongitude, WLlatitude, z3, 130, WLvelocity,'d', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Tide Gauge','LineWidth',0.1,'MarkerEdgeAlpha',0.6)
pbaspect([1 1 1])
legend('Location','southeast')
xlim([-125 -121.5])
ylim([45.5 49])
zlim([-3 6])
h = colorbar;
caxis([-3 3])
ylabel(h, 'Vertical Velocity (mm/yr)')
load PNWcoast.dat;
zz = zeros(1,length(PNWcoast));
zz = zz+6;
plot3(PNWcoast(:,1),PNWcoast(:,2), zz,'k', 'LineWidth', 1, 'DisplayName','Coastline')
plot3(contour20(:,1),contour20(:,2), zz20,'k', 'LineWidth', 2, 'DisplayName','20 km contour','LineStyle','-')
plot3(contour40(:,1),contour40(:,2), zz40,'k', 'LineWidth', 2, 'DisplayName','40 km contour','LineStyle','--')
plot3(contour60(:,1),contour60(:,2), zz60,'k', 'LineWidth', 2, 'DisplayName','60 km contour','LineStyle',':')
plot3(contour80(:,1),contour80(:,2), zz80,'k', 'LineWidth', 2, 'DisplayName','80 km contour','LineStyle','-.')
plot3(contour100(:,1),contour100(:,2), zz100,'k', 'LineWidth', 2, 'DisplayName','100 km contour','LineStyle','-')

% data map
figure
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('Data Used for Interpolation')
%scatter3(MODlongitude, MODlatitude, z4, 40, MODvelocity,'o', 'filled','DisplayName','Model')
scatter3(GPSlongitude, GPSlatitude, z1, 70, GPSvelocity,'<', 'filled','DisplayName','GPS')
scatter3(LEVlongitude, LEVlatitude, z2, 70, LEVvelocity,'s', 'filled','DisplayName','Leveling')
scatter3(WLlongitude, WLlatitude, z3, 130, WLvelocity,'d', 'filled','DisplayName','WL')
pbaspect([1 1 1])
legend('Location','southeast')
%xlim([-125 -121.5])
%ylim([45.5 49])
xlim([-125.5 -121])
ylim([45 50])
zlim([-3 6])
h = colorbar;
caxis([-3 3])
ylabel(h, 'Vertical Velocity (mm/yr)')
load PNWcoast.dat;
zz = zeros(1,length(PNWcoast));
zz = zz+6;
plot3(PNWcoast(:,1),PNWcoast(:,2), zz,'k', 'LineWidth', 2, 'DisplayName','Coastline')


% data density figure
qq = [];
for i = 1:46; %rows
    for ii = 1:51 %columns
        % find data that exists in this cell
        dref = find( (latlonvel(:,2) >= xq(i,ii)) & (latlonvel(:,2) <= (xq(i,ii)+0.09999999999)) & (latlonvel(:,1) >= yq(i,ii)) & (latlonvel(:,1) <= (yq(i,ii)+0.09999999999)) );
        % find marker type for this cell
        % 0 = no data in this cell
        % 1 = 1 data point in this cell
        % 2 = 2 or more data points of a single type in this cell
        % 3 = 2 or more data points of varying types in this cell
        dreflength = length(dref);
        dreftype = type6(dref);
        if dreflength >= 2;
            dreftypcmp = strcmp(dreftype(1,1), dreftype(:,1));
            dreftyplog = all(dreftypcmp(:) == 1);
        else
        end
        if dreflength == 0;
            qq(i,ii) = 0;
        elseif dreflength == 1;
            qq(i,ii) = 1;
        elseif dreflength >= 2 && dreftyplog == 1;
            qq(i,ii) = 2;
        else dreflength >= 2 && dreftyplog == 0;
            qq(i,ii) = 3;
        end
        dref = [];
    end
end

figure
surf(xq, yq, qq, 'EdgeColor','none') % interpolated vertical velocity surface
map = [1,0,0;1,0.5,0;1,1,0;0,1,0];
colormap(map)
xlabel('Longitude')
ylabel('Latitude')
zlabel('Data Indicator')
hold on
title('Data Quality Map')
xlim([-125 -121.5])
ylim([45.5 49])
zlim([0 6])
h = colorbar;
caxis([0 3])
ylabel(h, 'Data Quality')
set(h,'YTick',[0.375,1.125,1.875,2.625])
Labels = {'No Data','Single Data Point','>=2 Data Points, 1 Type','>=2 Data Types'};
set(h, 'TickLabels', Labels)
h.Location = 'westoutside';
pbaspect([1 1 1])
zz = zeros(1,length(PNWcoast));
zz = zz+6;
plot3(PNWcoast(:,1),PNWcoast(:,2), zz,'k', 'LineWidth', 2, 'DisplayName','Coastline')


%% E-W profiles
% data map for profile selection
close all
figure
coast = plot(PNWcoast(:,1),PNWcoast(:,2),'k', 'LineWidth', 2, 'DisplayName','Coastline')
hold on
scatter(longituderef, latituderef, 40, velocityref, 'o', 'filled');
colormap jet
caxis([-3 3])
colorbar
pbaspect([1 1 1])
xlim([-125 -121.5])
ylim([45.5 50])
[pind,xs,ys] = selectdata('selectionmode','lasso','ignore',coast) 

% profile generation for selection
profilelat = latituderef(pind);
profilelon = longituderef(pind);
profiletype = typeref(pind);
profilevel = velocityref(pind);
profileuncert = uncertref(pind);
index = find(contains(profiletype,'GPS'));
profileGPSlatitude = profilelat(index);
profileGPSlongitude = profilelon(index);
profileGPSvelocity = profilevel(index);
profileGPSuncert = profileuncert(index);
index = find(contains(profiletype,'Leveling'));
profileLEVlatitude = profilelat(index);
profileLEVlongitude = profilelon(index);
profileLEVvelocity = profilevel(index);
profileLEVuncert = profileuncert(index);
index = find(contains(profiletype,'WL'));
profileWLlatitude = profilelat(index);
profileWLlongitude = profilelon(index);
profileWLvelocity = profilevel(index);
profileWLuncert = profileuncert(index);
%index = find(contains(profiletype,'Model'));
%profileMODlatitude = profilelat(index);
%profileMODlongitude = profilelon(index);
%profileMODvelocity = profilevel(index);
figure
colormap jet
xlabel('Longitude')
ylabel('Velocity')
hold on
title('46.2 degree Profile')
xmin = min(profilelon);
xminorig = xmin;
xmin = round(xmin,1) - 0.1;
xmax = max(profilelon);
xmaxorig = xmax;
xmax = round(xmax,1) + 0.1;
xlim([xmin xmax])
set(gca,'XTick',xmin:0.1:xmax)
ymin = min(profilevel);
ymin = round(ymin,1) - 0.1;
ymax = max(profilevel);
ymax = round(ymax,1) + 0.1;
%ylim([ymin ymax])
ylim([-4 4])
legend('Location','southeast')
% simplified interpolated velocity for the selected points, assuming a straight line between the max and min longitude
x2 = xmaxorig;
x1 = xminorig;
xref = [xmin:0.1:xmax];
y2ref = find(profilelon == x2);
y2ref(2:100,1) = 1;
y2ref(2:100,:) = [];
y2lat = profilelat(y2ref);
y1ref = find(profilelon == x1);
y1ref(2:100,1) = 1;
y1ref(2:100,:) = [];
y1lat = profilelat(y1ref);
m = (y2lat-y1lat)/(x2-x1);
b = y2lat-(m*x2);
ylatref = (m.*xref)+b;
%plot(xref,F(xref,ylatref), 'LineWidth', 5, 'DisplayName','Interpolated Velocity')
%scatter(profileMODlongitude, profileMODvelocity, 40, 'o', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Model')
scatter(profileGPSlongitude, profileGPSvelocity, 70, '<', 'filled', 'MarkerEdgeColor', 'k','DisplayName','GPS')
scatter(profileLEVlongitude, profileLEVvelocity, 70, 's', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Leveling')
scatter(profileWLlongitude, profileWLvelocity, 130, 'd', 'filled', 'MarkerEdgeColor', 'k','DisplayName','WL')
hh = errorbar(profileGPSlongitude, profileGPSvelocity, profileGPSuncert, 'k','HandleVisibility','off');
hh.LineStyle = 'none';
jj = errorbar(profileLEVlongitude, profileLEVvelocity, profileLEVuncert, 'k','HandleVisibility','off');
jj.LineStyle = 'none';
kk = errorbar(profileWLlongitude, profileWLvelocity, profileWLuncert, 'k','HandleVisibility','off');
kk.LineStyle = 'none';
%% N-S profiles
% data map for profile selection
figure
coast = plot(washington_coastline(:,1),washington_coastline(:,2),'k', 'LineWidth', 2, 'DisplayName','Coastline')
hold on
scatter(longituderef, latituderef, 40, 'o', 'filled');
pbaspect([1 1 1])
xlim([-125 -121.5])
ylim([45.5 49])
[pind,xs,ys] = selectdata('selectionmode','lasso','ignore',coast) 

% profile generation for selection
profilelat = latituderef(pind);
profilelon = longituderef(pind);
profiletype = typeref(pind);
profilevel = velocityref(pind);
index = find(contains(profiletype,'GPS'));
profileGPSlatitude = profilelat(index);
profileGPSlongitude = profilelon(index);
profileGPSvelocity = profilevel(index);
index = find(contains(profiletype,'Leveling'));
profileLEVlatitude = profilelat(index);
profileLEVlongitude = profilelon(index);
profileLEVvelocity = profilevel(index);
index = find(contains(profiletype,'WL'));
profileWLlatitude = profilelat(index);
profileWLlongitude = profilelon(index);
profileWLvelocity = profilevel(index);
index = find(contains(profiletype,'Model'));
profileMODlatitude = profilelat(index);
profileMODlongitude = profilelon(index);
profileMODvelocity = profilevel(index);
figure
colormap jet
xlabel('Latitude')
ylabel('Velocity')
hold on
title('Profile for selected data')
xmin = min(profilelat);
xminorig = xmin;
xmin = round(xmin,1) - 0.1;
xmax = max(profilelat);
xmaxorig = xmax;
xmax = round(xmax,1) + 0.1;
xlim([xmin xmax])
set(gca,'XTick',xmin:0.1:xmax)
ymin = min(profilevel);
ymin = round(ymin,1) - 0.1;
ymax = max(profilevel);
ymax = round(ymax,1) + 0.1;
ylim([ymin ymax])
legend('Location','southeast')
% simplified interpolated velocity for the selected points, assuming a straight line between the max and min longitude
y2 = xmaxorig;
y1 = xminorig;
yref = [xmin:0.1:xmax];
x2ref = find(profilelat == y2);
x2ref(2:100,1) = 1;
x2ref(2:100,:) = [];
x2lon = profilelon(x2ref);
x1ref = find(profilelat == y1);
x1ref(2:100,1) = 1;
x1ref(2:100,:) = [];
x1lon = profilelon(x1ref);
m = (y2-y1)/(x2lon-x1lon);
b = y2-(m*x2lon);
xlonref = (yref-b)./m;
plot(yref,F(xlonref,yref), 'LineWidth', 5, 'DisplayName','Interpolated Velocity')
scatter(profileMODlatitude, profileMODvelocity, 40, 'o', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Model')
scatter(profileGPSlatitude, profileGPSvelocity, 70, '<', 'filled', 'MarkerEdgeColor', 'k','DisplayName','GPS')
scatter(profileLEVlatitude, profileLEVvelocity, 70, 's', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Leveling')
scatter(profileWLlatitude, profileWLvelocity, 130, 'd', 'filled', 'MarkerEdgeColor', 'k','DisplayName','WL')

%%

% generates data for coastline map

load('washington_coastline.mat');
latitude = washington_coastline(:, 2); 
longitude = washington_coastline(:, 1);
velocity = F(longitude,latitude);

figure
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('Interpolated Coastline Velocity with Model')
scatter(longitude, latitude,20, velocity,'.')
xlim([-125 -122])
ylim([46 49.3])
h = colorbar;
caxis([-3 3])
ylabel(h, 'Vertical Velocity (mm/yr)')
pbaspect([1 1 1])

% no model

E = scatteredInterpolant(longitude6, latitude6, velocity6);
velocity3 = E(longitude,latitude);
figure
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('Interpolated Coastline Velocity without Model')
scatter(longitude, latitude,20, velocity3,'.')
xlim([-125 -122])
ylim([46 49.3])
h = colorbar;
caxis([-3 3])
ylabel(h, 'Vertical Velocity (mm/yr)')
pbaspect([1 1 1])

[xq,yq] = meshgrid(-126:0.1:-121, 45:0.1:49.5);
vq = E(xq, yq);
figure
surf(xq, yq, vq, 'EdgeColor','none') % interpolated vertical velocity surface
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('Interpolated Velocity - no model')
xlim([-125 -121.5])
ylim([46 49])
zlim([-3 6])
h = colorbar;
caxis([-3 3])
ylabel(h, 'Vertical Velocity (mm/yr)')
pbaspect([1 1 1])
zz = zeros(1,length(washington_coastline));
zz = zz+6;
plot3(washington_coastline(:,1),washington_coastline(:,2), zz,'k', 'LineWidth', 2, 'DisplayName','Coastline')


% difference figure 
velocity4 = velocity3-velocity; % -1.5 to 2
figure
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('No Model - Model Velocity Difference')
scatter(longitude, latitude,20, velocity4,'.')
xlim([-125 -122])
ylim([46 49.3])
h = colorbar;
caxis([-1 1])
ylabel(h, 'Vertical Velocity (mm/yr)')
pbaspect([1 1 1])

%% figure for cumulative vertical displacement in 2100, relative to 2017
cumu83 = velocity.*83;
cumu83CM = cumu83./10;
cumu83M = cumu83CM./100;
figure
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Vertical Displacement')
hold on
title('Cumulative Vertical Displacement by 2100 (m)')
scatter(longitude, latitude,20, cumu83M,'.')
xlim([-125 -122])
ylim([46 49.3])
h = colorbar;
caxis([-1.0 1.0])
pbaspect([1 1 1])
ylabel(h, 'Vertical Displacement (m)')

% figure for baseline gradient subtraction (GIA)
%-123 to -122 averaging: vertical strip, mean taken horizontally, then
%subtracted from model surface of same latitude
[num] = xlsread('WAmodelpointsRay6.xlsx'); % reads data file and imports all elements that are numbers
latitudeLOC = num(:, 2); 
longitudeLOC = num(:, 1);
velocityLOC = num(:, 3);

F = scatteredInterpolant(longitudeLOC, latitudeLOC, velocityLOC, 'linear');
[xq,yq] = meshgrid(-126:0.1:-121, 45:0.1:49.5); % creates query points for interpolated plot
vq = F(xq, yq);


wq = vq;
vqsize = size(vq,1);
cattyTemp = [];
for i = 1:vqsize
    rowmean = mean(vq(i,31:41));
    wq(i,:) = wq(i,:)-rowmean;
    cattyTemp(i,1) = yq(i,1);
    cattyTemp(i,2) = rowmean;
end
fig = figure
surf(xq, yq, wq, 'EdgeColor','none') % interpolated vertical velocity surface with baseline gradient subtracted
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('Tectonic Model minus GIA')
xlim([-125 -121.5])
ylim([45.5 49])
zlim([-2 4])
h = colorbar;
caxis([-2 4])
%pbaspect([1 1 1])
ylabel(h, 'Vertical Velocity (mm/yr)')
print(fig, '-dpdf', '-r3000')

x2q = xq.';
y2q = yq.';
w2q = wq.';
G = griddedInterpolant(x2q,y2q,w2q)
noGiaVelocity = G(longitude,latitude);
% cumulative vertical displacement from 1700 to 2100, was previously 400 
cumu400 = noGiaVelocity.*525;
megathrust = cumu83-cumu400;
megathrustCM = megathrust./10;
megathrustM = megathrustCM./100
fig = figure
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Vertical Displacement')
hold on
title('Vertical displacement if rupture before 2100')
scatter(longitude, latitude,20, megathrustM,'.')
xlim([-125 -122])
ylim([46 49.3])
h = colorbar;
caxis([-1.0 1.0])
pbaspect([1 1 1])
ylabel(h, 'Vertical Displacement (m)')
%print(fig, '-dpdf', '-r3000')
%% generates data for  profiles

[num] = xlsread('profiledata.xlsx');
latitude = num(:, 1); 
longitude = num(:, 2);

interpvelarray = [];
for i=1:length(longitude)
[rowx,colx] = find(abs(abs(xq-longitude(i)) == min(abs(xq(:)-longitude(i)))));
[rowy,coly] = find(abs(abs(yq-latitude(i)) == min(abs(yq(:)-latitude(i)))));
interpVel = vq(rowy(1),colx(1));
interpvelarray = vertcat(interpvelarray, interpVel);
end

%% estimate of uncertainty
close all
clear all
format long
tic

[num, txt] = xlsread('new.xlsx'); % reads data file and imports all elements that are numbers
numNaN = any(isnan(num), 1); % finds NaN values in 'num'
num(:, numNaN) = []; % deletes columns with NaN values in 'num'

% exclude data of quality index 2 and 4
exclude2 = find(num(:,5) == 2);
num(exclude2,:) = [];
txt(exclude2,:) = [];
exclude4 = find(num(:,5) == 4);
num(exclude4,:) = [];
txt(exclude4,:) = [];

latitude = num(:, 2); 
longitude = num(:, 1);
velocity = num(:, 3);
standardDeviation = num(:, 4);
txt(1,:) = [];
station = txt(:, 1); 
type = txt(:, 6);

% first insert uncertainty in 4th column of modeldata.xlsx file
[num, txt] = xlsread('modeldata.xlsx');
modellatitude = num(:, 2); 
modellongitude = num(:, 1);
modelvelocity = num(:, 3);
modelstandarddeviation = num(:,4);
modeltype = txt(:,1);
type = vertcat(type, modeltype); 
latitude = vertcat(latitude, modellatitude); 
longitude = vertcat(longitude, modellongitude); 
velocity = vertcat(velocity, modelvelocity); 
standardDeviation = vertcat(standardDeviation, modelstandarddeviation); 
latitude2 = latitude;
longitude2 = longitude;
velocity2 = velocity;
standardDeviation2 = standardDeviation;

F = scatteredInterpolant(longitude, latitude, velocity);
[xq,yq] = meshgrid(-126:0.1:-121, 45:0.1:49.5); % creates query points for interpolated plot
points = [xq(:), yq(:)];
vq = F(xq, yq);


longitudeFil = longitude;
latitudeFil = latitude;
velocityFil = velocity;
E = scatteredInterpolant(longitudeFil, latitudeFil, velocityFil);
standardDeviationFil = standardDeviation;
latitude = latitude2;
longitude = longitude2;
velocity = velocity2;
standardDeviation = standardDeviation2;

rndStDev = [];
load('washington_coastline.mat');
%coastlatitude = points(:, 2); 
%coastlongitude = points(:, 1);
coastlatitude = washington_coastline(:, 2); 
coastlongitude = washington_coastline(:, 1);
% load('wa_coastline_tenth.mat');
% coastlatitude = wa_coastline_tenth(:, 2); 
% coastlongitude = wa_coastline_tenth(:, 1);

% loop through 1 million iterations
for i = 1:10000
	% generates random numbers from normal distributions (with mean specified by the 
	% values stored in *velocity* and the standard deviation specified by the values in *standardDeviation*)
    rndStDev(:,i) = normrnd(velocity,standardDeviation);
	% calculate a new height field based on the new randomly sampled vertical velocities
    F = scatteredInterpolant(longitude, latitude, rndStDev(:,i));
    vq = F(xq, yq);
        interpvelarray = [];
		% for each point on the coast:
        for ii=1:length(coastlongitude)
			% find the row and column of longitude grid closest POI on coast
            [rowx,colx] = find(abs(abs(xq-coastlongitude(ii)) == min(abs(xq(:)-coastlongitude(ii)))));
			% find the row and column of latitude grid closest POI on coast
            [rowy,coly] = find(abs(abs(yq-coastlatitude(ii)) == min(abs(yq(:)-coastlatitude(ii)))));
            % extract velocity from interpolation at coast location
			interpVel = vq(rowy(1),colx(1));
            interpvelarray = vertcat(interpvelarray, interpVel);
        end
    coastCat(:,i) = interpvelarray;
    i
end

coastUncertRange = range(coastCat,2);      
coastUncertMean = mean(coastCat,2);        
coastUncertStd = std(coastCat,0,2);

saveme(:,1) = coastlongitude;
saveme(:,2) = coastlatitude;
saveme(:,3) = E(coastlongitude,coastlatitude);
saveme(:,4) = coastUncertStd;
%dlmwrite('coast_data_newCoastFile_NEW.txt', saveme, 'delimiter', ' ','precision',7)
dlmwrite('interp_uncert_Water.txt', saveme, 'delimiter', ' ','precision',7)

U = scatteredInterpolant(coastlongitude, coastlatitude, coastUncertStd);
uq = U(xq, yq);
uq = abs(uq);
arcgridwrite('arc_uncert_surf_Water.asc',xq(1,:),yq(:,1),uq);

toc
%% uncert fig
%saveme = load('coast_data_newCoastFile_NEW.txt');
saveme = load('interp_uncert_Water.txt');
coastlongitude = saveme(:,1);
coastlatitude = saveme(:,2);
E = saveme(:,3);
coastUncertStd = saveme(:,4);

figure
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Velocity')
hold on
title('Interpolated Coastline Velocity with Model')
scatter(coastlongitude, coastlatitude,20, E,'.')
xlim([-125 -122])
ylim([46 49.3])
h = colorbar;
caxis([-2 4])
ylabel(h, 'Vertical Velocity (mm/yr)')
pbaspect([1 1 1])

% uncertainty figures

figure
colormap jet
xlabel('Longitude')
ylabel('Latitude')
zlabel('Uncertainty Standard Deviation')
hold on
title('Interpolated Coastline Velocity Standard Deviation')
scatter(coastlongitude, coastlatitude, 20, coastUncertStd,'.')
xlim([-125 -122])
ylim([46 49.3])
h = colorbar;
caxis([0 4.0])
pbaspect([1 1 1])
ylabel(h, 'Vertical Velocity Standard Deviation (mm/yr)')

% standardDeviation = abs(standardDeviation);
% index = find(contains(type,'GPS'));
% GPSlatitude = latitude(index);
% GPSlongitude = longitude(index);
% GPSstd = standardDeviation(index);
% index = find(contains(type,'Leveling'));
% LEVlatitude = latitude(index);
% LEVlongitude = longitude(index);
% LEVstd = standardDeviation(index);
% index = find(contains(type,'WL'));
% WLlatitude = latitude(index);
% WLlongitude = longitude(index);
% WLstd = standardDeviation(index);
% index = find(contains(type,'Model'));
% MODlatitude = latitude(index);
% MODlongitude = longitude(index);
% MODstd = standardDeviation(index);
% z = zeros(length(latitude),1);
% z = z+8.5;
% z1 = zeros(length(GPSlatitude),1);
% z1 = z1+8.5;
% z2 = zeros(length(LEVlatitude),1);
% z2 = z2+8.5;
% z3 = zeros(length(WLlatitude),1);
% z3 = z3+8.5;
% z4 = zeros(length(MODlatitude),1);
% z4 = z4+8.5;
% E = scatteredInterpolant(longitude, latitude, standardDeviation);
% [xxq,yyq] = meshgrid(-126:0.1:-121, 45:0.1:49.5); % creates query points for interpolated plot
% vvq = E(xxq, yyq);
% vvq = abs(vvq)
% figure
% surf(xxq, yyq, vvq, 'EdgeColor','none','DisplayName', 'Interpolated Uncertainty')
% colormap jet
% xlabel('Longitude')
% ylabel('Latitude')
% zlabel('Standard Deviation')
% hold on
% title('Interpolated Standard Deviation')
% scatter3(GPSlongitude, GPSlatitude, z1, 40, GPSstd,'o', 'filled', 'MarkerEdgeColor', 'k','DisplayName','GPS')
% scatter3(LEVlongitude, LEVlatitude, z2, 40, LEVstd,'s', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Leveling')
% scatter3(WLlongitude, WLlatitude, z3, 40, WLstd,'d', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Tide Gauge')
% %scatter3(MODlongitude, MODlatitude, z4, 40, MODstd,'<', 'filled', 'MarkerEdgeColor', 'k','DisplayName','Model')
% pbaspect([1 1 1])
% legend('Location','southeast')
% xlim([-125 -121.5])
% ylim([45.5 49])
% zlim([0 9])
% h = colorbar;
% caxis([0 4.0])
% pbaspect([1 1 1])
% ylabel(h, 'Vertical Velocity (mm/yr)')
% load PNWcoast.dat;
% zz = zeros(1,length(PNWcoast));
% zz = zz+8.7;
% plot3(PNWcoast(:,1),PNWcoast(:,2), zz,'k', 'LineWidth', 2)
% % writes .asc file
% %arcgridwrite('arc_uncert_surf_new.asc',xxq(1,:),yyq(:,1),vvq);
% 
% % Points of Interest histograms
% [num] = xlsread('POI.xlsx');
% POI(:,1) = num(:, 1); %lat
% POI(:,2) = num(:, 2); %long
% REF(:,1) = coastlatitude;
% REF(:,2) = coastlongitude;
% IDX = knnsearch(REF,POI);
% POIpoints = coastCat(IDX,:);
% 
% figure
% neahbay = histogram(POIpoints(1,:),'BinWidth',0.2)
% xlabel('Vertical Velocity (mm/yr)')
% ylabel('Number of Occurrences')
% title('Neah Bay Uncertainty Histogram - With Model')
% xlim([-10 10])
% ylim([0 200])
% 
% figure
% portangeles = histogram(POIpoints(2,:),'BinWidth',0.2)
% xlabel('Vertical Velocity (mm/yr)')
% ylabel('Number of Occurrences')
% title('Port Angeles Uncertainty Histogram - With Model')
% xlim([-10 10])
% ylim([0 200])
% 
% figure
% seattle = histogram(POIpoints(3,:),'BinWidth',0.2)
% xlabel('Vertical Velocity (mm/yr)')
% ylabel('Number of Occurrences')
% title('Seattle Uncertainty Histogram - With Model')
% xlim([-10 10])
% ylim([0 200])
% 
% figure
% tokepoint = histogram(POIpoints(4,:),'BinWidth',0.2)
% xlabel('Vertical Velocity (mm/yr)')
% ylabel('Number of Occurrences')
% title('Toke Point Uncertainty Histogram - With Model')
% xlim([-10 10])
% ylim([0 200])
% 
% figure
% cherrypoint = histogram(POIpoints(5,:),'BinWidth',0.2)
% xlabel('Vertical Velocity (mm/yr)')
% ylabel('Number of Occurrences')
% title('Cherry Point Uncertainty Histogram - With Model')
% xlim([-10 10])
% ylim([0 200])
% 
% figure
% portrenfew = histogram(POIpoints(6,:),'BinWidth',0.2)
% xlabel('Vertical Velocity (mm/yr)')
% ylabel('Number of Occurrences')
% title('Port Renfrew Uncertainty Histogram - With Model')
% xlim([-10 10])
% ylim([0 200])
% 
% figure
% pointofinterest = histogram(POIpoints(7,:),'BinWidth',0.2)
% xlabel('Vertical Velocity (mm/yr)')
% ylabel('Number of Occurrences')
% title('47.8, -124.5 Uncertainty Histogram - With Model')
% xlim([-10 10])
% ylim([0 200])

%% script to round coastline file
load('washington_coastline.mat');
wa_coastline_tenth = round(washington_coastline, 1);
wa_coastline_tenth = unique(wa_coastline_tenth, 'rows');
