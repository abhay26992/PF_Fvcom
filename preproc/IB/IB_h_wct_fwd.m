% Read the Ice-Bridge bathymetry data for Petermann

%%
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;


%% Ice - Bridge coordinates:

ice_bridge=readtable('petermann_tinto_inv_bathy.csv');

ib_lon=table2array(ice_bridge(:,2));ib_lat=table2array(ice_bridge(:,3));

a = 6378137;
e = 0.08181919;
phi_c = 70;
lambda_0 = -45;
[ib_x,ib_y]=polarstereo_fwd(ib_lat,ib_lon,a,e,phi_c,lambda_0);

ib_line=table2array(ice_bridge(:,1));ib_bth=table2array(ice_bridge(:,4)).*(-1);

figure(1);clf
scatter(ib_x,ib_y,40,ib_bth);colorbar;colormap jet; %or use lon-lat
load ./M
hold on;
scatter(Mobj.x,Mobj.y,1,'k');

load Lon_Crack;load Lat_Crack;

figure(2);clf
plot_field(Mobj,Mobj.h,'coordinate','spherical');colormap jet;colorbar;
caxis([0 1100]);
hold on;scatter(ib_lon,ib_lat,40,ib_bth);
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
draft_ibg = interp2(BM_x,BM_y,pm_draft,ib_x,ib_y);
wct_ibg=ib_bth-draft_ibg;

figure(3);clf
plot_field(Mobj,Mobj.h-Mobj.zisf,'coordinate','spherical');colormap jet;colorbar;
caxis([0 900]);
hold on;scatter(ib_lon,ib_lat,80,wct_ibg);
scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');