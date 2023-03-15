%% Generate the new cracked geometry for FVCOM 

% Steps:

% 1) Build virtual rasters using the R(04), G(03), B(02) bands for each tile 
% from the 10-m Sentinel-2 data. Merge the raster tiles.

% 2) Following the Petermann crack on the tile, create a 'crack' shapefile 
% and extract the coordinates. 

% 3) Matlab --> Read the coordinates of the crack (x,y in UTM Zone 21 N) 
% and convert to the FVCOM spherical coordinates.

% 4) Matlab --> Build a polygon using the "crack coordinates" around the 
% FVCOM Petermann ice-shelf. 

% 5) Set Zisf=0 in this polygon region and write the new geometry.

%% Step 3:

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;

Crack=csvread('Crack_XY.csv',1,0);
x=Crack(:,1);y=Crack(:,2);

UTM_Zone='21 N';
UTM_Zone=repmat(UTM_Zone,length(x),1);
[Lat,Lon]=utm2deg(x,y,UTM_Zone);

%% Step 4:


%Load the inner-sill Mobj
load ./IBN_M.mat

% Visualise the crack over 'zisf'
figure(1);clf
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);
hold on;
scatter(Lon,Lat,20,'k','o','filled');

% Account for Sentinel-2 (2019) data vs FVCOM coastline/domain mismatch to
% complete the crack:

% 1. Extend the points on the right flank 
Lon_ext=[-60.77;-60.76;-60.74];
Lat_ext=[80.98;80.98;80.98];

scatter(Lon_ext,Lat_ext,20,'r','o','filled');

%2. Append these points to the "Crack":

% Note: The order is important. The crack coordinates (Lon,Lat) start from
% the right flank and end at the left flank. Therefore, points added on the
% edge of the right flank are the first Lon,Lat pairs and must be added
% before the Lon,Lat array starts.

Lon_Crack=[Lon_ext;Lon];
Lat_Crack=[Lat_ext;Lat];

% Do a quick-save:
% save Lon_Crack Lon_Crack; save Lat_Crack Lat_Crack;

% Visualise again:
figure(2);clf
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);
hold on;
scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');


%%%%%%%%%%%%%%%%%%%%%% Close the 'Crack' Polygon %%%%%%%%%%%%%%%%%%%%%%%%%%

xclose=[-60.75;-61.5;-62.04;-62.5;-62;-61.6];
yclose=[81;81.25;81.26;81.15;81;80.83];

% Append the indices to the Crack Coordinates:

% Note: The order is important. We start with xclose(1),yclose(1) which 
% and in a counter-clockwise manner approach the left Petermann flank. We, 
% therefore, need to append our Crack coordinates [Lon_Crack,Lat_Crack]
% after these points to close the polygon perfectly.

xpol=[xclose;flipud(Lon_Crack)]; %flipud: order -> left flank to right
ypol=[yclose;flipud(Lat_Crack)];

% Visualise:
figure(3);clf
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);
hold on;
scatter(xpol,ypol,20,'k','o','filled');

indx =find(inpolygon(Mobj.lon,Mobj.lat,xpol,ypol) & Mobj.zisf>0);
Mobj.zisf(indx)=0;

% Visualise:
figure(4);clf
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);
hold on;scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');

%save ./crackIBNM.mat Mobj %Save the Mobj with the new cracked "zisf"

%% Smoothing:

% 1. Extend the Crack Northward:
Ext_Lat_Crack=Lat_Crack+0.032;

% Visualise:
figure(5);clf
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);hold on;
scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');
scatter(Lon_Crack,Ext_Lat_Crack,20,'r','o','filled');


% 2. Find Mobj.zisf values at Lon_Crack and Lat_Crack:

% At the crack nodes (variable = ind), some are set to 0 and some have
% values. As all the nodes nortward of the crack are set to 0, we get steep
% gradients. Therefore, first we set the value at the crack line to the
% original "zisf" values, and finally, smooth it to a desired minimum value
% at the front, slightly northward of the crack line.

clear Mobj; %clear the Mobj with the cracked geometry
load ./IBN_M.mat; %load the original Mobj
ind=NaN(size(Lon_Crack));

for i=1:length(Lon_Crack)
    l=spheredist(Mobj.lat,Mobj.lon,Lat_Crack(i,1),Lon_Crack(i,1));
    ind(i,1)=find(l==min(l));
end

fv_lon=Mobj.lon(ind);fv_lat=Mobj.lat(ind);
fv_ice=Mobj.zisf(ind);

figure(6);clf
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);hold on;
scatter(fv_lon,fv_lat,100,fv_ice);


% 3. Interpolate:

%%%%%%%%%%%%% a. Allocate fv_ice values at the crack %%%%%%%%%%%%%%%%%%%%%%

org_Mobj=Mobj;clear Mobj; %Original geometry -> org_Mobj
load ./crackIBNM.mat; % Load the crack geometry -> Mobj


figure(7);clf      % Crack before allocation
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);hold on;


Mobj.zisf(ind)=org_Mobj.zisf(ind); % Allocate values

figure(8);clf      % Crack after allocation
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);hold on;


%%%%%%%%%%%%%%%%%%%%%%%% b. Get a Polygon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Following the comments from the xpol,ypol generation from above:

xp_r=[-60.7;-60.7];yp_r=[80.99;81];
xp_l=[-61.6;-61.6];yp_l=[80.83;80.84];

smo_polx=[Lon_Crack;xp_l;flipud(Lon_Crack);xp_r]; 
smo_poly=[Lat_Crack;yp_l;flipud(Ext_Lat_Crack);yp_r];

%Test the Polygon:
test_i =find(inpolygon(Mobj.lon,Mobj.lat,smo_polx,smo_poly) & Mobj.zisf==0);
var_temp=Mobj.zisf;
var_temp(test_i)=5000; %Assign an absurdly high dummy value within the polygon

%Visualise
figure(9);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 400]);hold on;
scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');
scatter(Lon_Crack,Ext_Lat_Crack,20,'r','o','filled');


%%%%%%%%%%%%%%%%%%%%%%% c. Provide a smooth-ramp %%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the lat,(lon) and zisf values at the external crack node points:
ext_ind=NaN(size(Lon_Crack));

for i=1:length(Lon_Crack)
    l=spheredist(Mobj.lat,Mobj.lon,Ext_Lat_Crack(i,1),Lon_Crack(i,1));
    ext_ind(i,1)=find(l==min(l));
end

fv_ext_lat=Mobj.lat(ext_ind);fv_ext_ice=Mobj.zisf(ext_ind);
fv_ext_ice=fv_ext_ice+50; %Set the value at the external line to 50

%Find the indices of the node points that fall between the two crack lines:
ind_inside_lines=find(inpolygon(Mobj.lon,Mobj.lat,smo_polx,smo_poly));  
inside_lon=Mobj.lon(ind_inside_lines);inside_lat=Mobj.lat(ind_inside_lines);

figure(10);clf % Visualise:
scatter(fv_lon,fv_lat,50,fv_ice);hold on;
scatter(fv_lon,fv_ext_lat,50,fv_ext_ice);
colorbar;colormap jet;
scatter(inside_lon,inside_lat,1,'k','filled');

% Create a variable which contains the zisf values at the two bounding
% crack lines and use it to create the interpolant:

ice_bound=[fv_ice;0;0;fv_ext_ice;0;0]; % 0 -> ice-value (outside the domain) at the x-y points that help close the polygon
F=scatteredInterpolant(smo_polx,smo_poly,ice_bound);

%Obtain the values at the nodes inside the two bounding lines using the
%interpolant:
ice_val=F(inside_lon,inside_lat);

%Visualise:
figure(11);clf
scatter(inside_lon,inside_lat,10,ice_val);colorbar;colormap jet;
hold on;
scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');
scatter(Lon_Crack,Ext_Lat_Crack,20,'r','o','filled');

% Append the values to a new "zisf" variable:
zisf=Mobj.zisf;
zisf(ind_inside_lines)=ice_val;

%Visualise:
figure(12);clf 
plot_field(Mobj,zisf,'coordinate','spherical');
colormap jet; caxis([0 400]);
hold on;
scatter(inside_lon,inside_lat,10,'k','filled');

% Clean - up noise after interpolation:
%[noise_x,noise_y]=ginput;
%save noise_x noise_x; save noise_y noise_y;

load noise_x;load noise_y;

scatter(noise_x,noise_y,10,'r','filled');

i_noise =find(inpolygon(Mobj.lon,Mobj.lat,noise_x,noise_y) & zisf>100);
zisf(i_noise)=100;

%Visualise:
figure(13);clf 
plot_field(Mobj,zisf,'coordinate','spherical');
colormap jet; caxis([0 400]);
hold on;
scatter(Lon_Crack,Lat_Crack,20,'k','o','filled');
scatter(Lon_Crack,Ext_Lat_Crack,20,'r','o','filled');

% Append the values to Mobj.zisf and save as smo_crackM.mat:
Mobj.zisf=zisf; 
%save ./smo_crackIBNM.mat Mobj;

%% Smooth the WCT near the front after making the crack:

%Similar procedure as above:

%1. Make the North-South lines. Since there are no 'crack' coordinates
%available. All points need to be manually generated.

%2. Get WCT values at these two lines and generate the polygon. Test the 
%polygon for leaks (fill with zeros).  

%3. Use Scattered Interpolant to smoothly interpolate - from south line WCT
%value to north line WCT value. Note: We do not set the value of the north
%line to a scalar value, because this is the WCT, not the zisf, so we need
%the WCT values from BM dataset, just that it needs to be smoothed out.

%4. Note that all this does is gradually arrive at the deeper WCT at the
%north line.

% %%
% 
% %1. Make the lines: ginput
% 
% figure(14);clf
% plot_field(Mobj,Mobj.h-Mobj.zisf,'coordinate','spherical');caxis([0 800]);shading flat;colormap jet; %h-> north points. zisf -> south points.
% hold on;
% % [wct_nx,wct_ny]=ginput;
% % [wct_sx,wct_sy]=ginput;
% % save wct_nx wct_nx;save wct_ny wct_ny;
% % save wct_sx wct_sx;save wct_sy wct_sy;
% 
% % load ./wct_nx; load ./wct_ny;
% % load ./wct_sx; load ./wct_sy;
% scatter(wct_nx,wct_ny,'k','o','filled');
% scatter(wct_sx,wct_sy,'k','o','filled');
% 
% %2.
% %2a. Get the lat,(lon) and wct values at the northern node points:
% wct_n_ind=NaN(size(wct_nx));
% 
% for i=1:length(wct_nx)
%     l=spheredist(Mobj.lat,Mobj.lon,wct_ny(i,1),wct_nx(i,1));
%     wct_n_ind(i,1)=find(l==min(l));
% end
% 
% wct_n_lon=Mobj.lon(wct_n_ind);wct_n_lat=Mobj.lat(wct_n_ind);
% wct_n_ice=Mobj.zisf(wct_n_ind);wct_n_h=Mobj.h(wct_n_ind);
% wct_n_wct=wct_n_h-wct_n_ice;
% 
% 
% %2b. Get the lat,(lon) and wct values at the southern node points:
% wct_s_ind=NaN(size(wct_sx));
% 
% for i=1:length(wct_sx)
%     l=spheredist(Mobj.lat,Mobj.lon,wct_sy(i,1),wct_sx(i,1));
%     wct_s_ind(i,1)=find(l==min(l));
% end
% 
% wct_s_lon=Mobj.lon(wct_s_ind);wct_s_lat=Mobj.lat(wct_s_ind);
% wct_s_ice=Mobj.zisf(wct_s_ind);wct_s_h=Mobj.h(wct_s_ind);
% wct_s_wct=wct_s_h-wct_s_ice;
% 
% %2c. Left right padding - not needed
% % wct_xl=[-61.6;-61.6];wct_yl=[80.83;80.84];
% % wct_xr=[-60.7;-60.7];wct_yr=[80.99;81];
% 
% smo_wctpolx=[wct_sx;flipud(wct_nx)]; 
% smo_wctpoly=[wct_sy;flipud(wct_ny)];
% 
% %Test the Polygon:
% test_wct_i =find(inpolygon(Mobj.lon,Mobj.lat,smo_wctpolx,smo_wctpoly));
% var_wct_temp=Mobj.zisf;
% var_wct_temp(test_wct_i)=5000; %Assign an absurdly high dummy value within the polygon
% 
% %Visualise
% figure(15);clf
% plot_field(Mobj,var_wct_temp,'coordinate','spherical');colormap jet;caxis([0 400]);hold on;
% scatter(wct_sx,wct_sy,20,'k','o','filled');
% scatter(wct_nx,wct_ny,20,'r','o','filled');
% 
% 
% %3. Smooth interpolation
% 
% %3a. Find the indices of the node points that fall between the two crack lines:
% ind_wct_inside_lines=find(inpolygon(Mobj.lon,Mobj.lat,smo_wctpolx,smo_wctpoly));  
% inside_wct_lon=Mobj.lon(ind_wct_inside_lines);inside_wct_lat=Mobj.lat(ind_wct_inside_lines);
% 
% figure(16);clf % Visualise:
% scatter(wct_s_lon,wct_s_lat,50,wct_s_wct);hold on;
% scatter(wct_n_lon,wct_n_lat,50,wct_n_wct);
% colorbar;colormap jet;
% scatter(inside_wct_lon,inside_wct_lat,1,'k','filled');
% 
% %3b. Create a variable which contains the wct values at the two bounding 
% % crack lines and use it to create the interpolant:
% 
% wct_bound=[wct_s_wct;wct_n_wct]; 
% F_wct=scatteredInterpolant(smo_wctpolx,smo_wctpoly,wct_bound);
% 
% %3c. Obtain the values at the nodes inside the two bounding lines using the
% % interpolant:
% wct_val=F_wct(inside_wct_lon,inside_wct_lat);
% 
% %Visualise:
% figure(17);clf
% scatter(inside_wct_lon,inside_wct_lat,10,wct_val);colorbar;colormap jet;
% hold on;
% scatter(wct_sx,wct_sy,20,'k','o','filled');
% scatter(wct_nx,wct_ny,20,'r','o','filled');
% 
% %3d. Append the values from "wct_val" to a new "my_wct" variable:
% my_wct=Mobj.h-Mobj.zisf;
% my_wct(ind_wct_inside_lines)=wct_val;
% %Remember to correct the 'h' now that we have changed the 'wct':
% my_h=my_wct+Mobj.zisf; %my_h (and my_wct) are the new geometry + Mobj.zisf (from before)
% 
% %Visualise:
% figure(18);clf 
% plot_field(Mobj,my_h,'coordinate','spherical'); %also check my_h
% colormap jet; caxis([0 800]);
% hold on;
% scatter(inside_wct_lon,inside_wct_lat,10,'k','filled');
% 
% %3e. Save h (if new wct is created)
% % Mobj.h=my_h; 
% % save ./smo_crackIBNM_wct.mat Mobj;
% 
% 
% % Notes:
% % This does not seem necessary. If needed: Figure (18) suggests that it is
% % best to make 2 small polygons on either flank and smooth (locally), the
% % patches on the flanks that show abrupt change. Smoothing the whole
% % section bound in here:
% % smo_wctpolx=[wct_sx;flipud(wct_nx)];smo_wctpoly=[wct_sy;flipud(wct_ny)];
% % causes errors, and is also not required, as the middle sectors are okay.
% 
% % Another suggestion is to use the Lon/Lat Crack from the zisf portion,
% % i.e. (Lon_Crack,Lat_Crack; Lon_Crack,Ext_Lat_Crack) instead of using
% % ginput to create (wct_nx,wct_ny; wct_sx,wct_sy). Since wct needs to be
% % corrected in this patch, using the former makes sense.
% 
% 
%% Create new input files 

%Set up model specific components:
Mobj = setup_metrics(Mobj);

casename=get_cstr;

% Write for FVCOM:
% switch to sphercial coordinates
Mobj.nativeCoords='spherical';
% dump mesh and connectivity
write_FVCOM_grid(Mobj,'../input/2019_Crack/grd.dat');
% dump bathymetry
write_FVCOM_bath(Mobj,'../input/2019_Crack/dep.dat');
% % dump open boundary node list
write_FVCOM_obc(Mobj,'../input/2019_Crack/obc.dat');
% dump sponge layer file
write_FVCOM_sponge(Mobj,'../input/2019_Crack/spg.dat');
% dump Coriolis file
write_FVCOM_cor(Mobj,'../input/2019_Crack/cor.dat');

fileprefix='pf_o_1';
write_FVCOM_iceshelfdraft(Mobj, fileprefix);

%% Plot to check:

clear all;close all;clc;
load ./smo_crackIBNM.mat;

figure(5);clf
plot_field(Mobj,Mobj.zisf);caxis([0 300]);shading flat;colormap jet;

figure(6); clf
plot_field(Mobj,Mobj.h);caxis([0 800]);shading flat;colormap jet;

figure(7); clf
plot_field(Mobj,Mobj.h-Mobj.zisf);caxis([0 800]);shading flat;colormap jet;