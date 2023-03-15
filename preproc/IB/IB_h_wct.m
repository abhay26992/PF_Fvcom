% Read the Ice-Bridge bathymetry data for Petermann

%%
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;

%% Bathymetry data:

ice_bridge=csvread('IGBTH4_20140207.csv',10,0);

% Refer to https://nsidc.org/data/IGBTH4/versions/1 for file structure, or 
% read the csv document.Use the Petermann Glacier ID (193) from the Line ID
% (first column) to trim the dataset:

trim_ind=find(ice_bridge(:,1)>=193 & ice_bridge(:,1)<194); 

% From the Tinto paper:
%OIB has flown six survey lines along the axis of Petermann Glacier and 
%fjord (Fig.3a). Two lines flown on 7 May 2011 surveyed the eastern (A–A')
%and western (B–B') sides of the Petermann fjord and continued across open
%water(Fig.3). Three shorter lines, flown on 24th March and 20th April 2010
%,without a magnetometer, are clustered on the eastern side of the fjord 
%and record smaller scale spatial variability. The maximum separation of 
%these lines is approximately 2.5 km, with one repeating the track of line 
%A–A' to within 100 m lateral separation.

line_id=ice_bridge(trim_ind,1); % We only have one of the 2011 tracks

lon=ice_bridge(trim_ind,4);lat=ice_bridge(trim_ind,5);
x=ice_bridge(trim_ind,6);y=ice_bridge(trim_ind,7);
bth=ice_bridge(trim_ind,end);bth=bth.*(-1);

figure(1);clf
scatter(x,y,40,bth);colorbar;colormap jet; %or use lon-lat
load ./M
hold on;
scatter(Mobj.x,Mobj.y,1,'k');

load Lon_Crack;load Lat_Crack;

figure(2);clf
plot_field(Mobj,Mobj.h,'coordinate','spherical');colormap jet;colorbar;
caxis([0 1200]);
hold on;scatter(lon,lat,80,bth);
scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');

%% Plot the "WCT" for the Ice-Bridge line:

% Get the ice-draft values for the line:

bm_x=ncread('BedMachineGreenland-2017-09-20','x'); %Cartesian X (m)
bm_y=ncread('BedMachineGreenland-2017-09-20','y'); %Cartesian Y (m)
bm_x=double(bm_x);bm_y=double(bm_y);[BM_x,BM_y] = meshgrid(bm_x,bm_y); 

ice_elev=ncread('BedMachineGreenland-2017-09-20','surface'); %Ice-Surface elevation (m)
thickness=ncread('BedMachineGreenland-2017-09-20','thickness'); %Ice-thickness (m)
ice_elev=double(ice_elev)';thickness=double(thickness)';
ice=thickness-ice_elev; %Draft = Thickness - Elevation

ii  = 1:5000;
jj = 1:5000;
% Remove all the ice that is not Petermann
pmii =2300:3200;
pmjj =1800:3000;
pm_draft = ice*0;
pm_draft(pmjj,pmii)=ice(pmjj,pmii);
pm_draft(1800:2300,2700:3200)=0; %Further cropping of the adjacent fjord

% Get the ice-draft values over the Ice-Bridge line:
draft_ibg = interp2(BM_x,BM_y,pm_draft,x,y);
wct_ibg=bth-draft_ibg;

figure(3);clf
plot_field(Mobj,Mobj.h-Mobj.zisf,'coordinate','spherical');
hold on;scatter(lon,lat,80,wct_ibg);
c=colorbar;c.Label.String='WCT [m]';colormap jet;caxis([0 1200]); %800 if exaggeration is needed
xlabel('Longitude');ylabel('Latitude');
scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');
