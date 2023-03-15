%% Calculate Subglacial Discharge @ Petermann using RACMO

% The drainage basin mask sits on the BedMachine grid. The Petermann basin 
% is all pixels with a value of 78.  

% Script that converts the BM grid x,y to lon,lat and finds indices that 
% correspond to the Petermann basin. These are used to find the RACMO nodes
% that fall within the basin. 

% Abhay Prakash (abhay.prakash@natgeo.su.se)

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;
%% Load the drainage basins and crop out NW Greenland region

load drainagebasins.mat;load IBN_M.mat;

figure(1);clf;
imagesc(x,y,drainagebasins);axis xy equal;hold on;
scatter(Mobj.x,Mobj.y,10,'k');

figure(2);clf;
kk  = 1:6800;
mm = 1:6800;
imagesc(x(kk),y(mm),drainagebasins(mm,kk)); axis xy equal; colorbar; 

x=x(kk);y=y(mm);db=drainagebasins(mm,kk);

%% Get the PMF drainage basin lat,lon

[x,y]=meshgrid(x,y);
[lat,lon]=polarstereo_inv(x,y);

pmdb_ind=find(db==78);
pm_db=db(pmdb_ind);
pmdb_lat=lat(pmdb_ind);pmdb_lon=lon(pmdb_ind);


%% Load racmo data

% Replace the 2-m air temperature with runoff when it becomes available.
% Also, the time dimenson (3-hourly time steps) needs to be trimmed down
% from the 2010-2018 period to the summer-months of the year of choice.

cdir=pwd; cd E:/Studies/Forcing/Atmospheric_Forcing/RACMO/Datasets
fn='t2m.nc';
ncdisp(fn)
rclat=ncread(fn,'lat');rclon=ncread(fn,'lon');t2m=ncread(fn,'t2m');
cd(cdir);

t2m=mean(squeeze(t2m),3); %Taking time-mean for now


%% Create a polygon (boundary) of the Petermann catchment

% for i=1:length(pmdb_lat)
%     for j=1:size(rclat,2)
%         dist(:,j)=spheredist(pmdb_lat(i),pmdb_lon(i),rclat(:,j),rclon(:,j));
%     end
%     ind(i) = find(dist == minmin(dist));
% end
% 
% ind = unique(ind,'stable');


% Using spheredist will become demanding, get the polygon instead:
bry=boundary(double(pmdb_lon),double(pmdb_lat),1); %Shrink factor = 0 gives the convex hull, here using 1

% Visualise
figure(3);clf
scatter(pmdb_lon,pmdb_lat,1,'k');hold on;
plot(pmdb_lon(bry),pmdb_lat(bry),'r','LineWidth',2);

%% Get the RACMO lat,lon indices falling within the Petermann basin

in_bry=inpolygon(rclon,rclat,pmdb_lon(bry),pmdb_lat(bry));
rclon_db=rclon(in_bry);rclat_db=rclat(in_bry);t2m_db=t2m(in_bry);

% Visualise
figure(4);clf
scatter(pmdb_lon,pmdb_lat,1,'k');hold on;
plot(pmdb_lon(bry),pmdb_lat(bry),'r','LineWidth',2);
scatter(rclon_db,rclat_db,10,'b');  

figure(5);clf
pcolor(rclon,rclat,t2m);colormap jet;colorbar;hold on;
scatter(pmdb_lon,pmdb_lat,2,'m','filled');
plot(pmdb_lon(bry),pmdb_lat(bry),'k','LineWidth',3);
scatter(rclon_db,rclat_db,5,'y','filled');
scatter(Mobj.lon,Mobj.lat,3,'r','filled');

%% Generate the product

% .....