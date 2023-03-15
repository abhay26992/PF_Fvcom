%% Gen_Sill_Geom

% Generate the sill geometry for PF Topography Experiment:

%%
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;

load ./M.mat
load ./Lon_Crack;load ./Lat_Crack;

%% 1. Extend the crack Northward and Southward

Ext_Lat_Crack_N=Lat_Crack+0.055;Ext_Lon_Crack_N=Lon_Crack-0.1;
Ext_Lat_Crack_S=Lat_Crack-0.20;Ext_Lon_Crack_S=Lon_Crack+0.55;

%% Make the large N-S polygon and set all h inside polygon = 0

% See Gen_Crack_Geom.m for detailed comments

x_fill_l=[-61.2;-61.4;-61.6;-61.7];y_fill_l=[80.62;80.7;80.8;80.86];
x_fill_r=[-60.1;-60.4;-60.6];y_fill_r=[80.85;80.95;81];


xp=[Ext_Lon_Crack_S;x_fill_l;flipud(Ext_Lon_Crack_N);x_fill_r];
yp=[Ext_Lat_Crack_S;y_fill_l;flipud(Ext_Lat_Crack_N);y_fill_r];


%Test the Polygon:
test_i =find(inpolygon(Mobj.lon,Mobj.lat,xp,yp));
var_temp=Mobj.h;
var_temp(test_i)=5000; %Test: Assign an absurdly high dummy value within the polygon

%Visualise
figure(2);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 1200]);hold on;

clear var_temp
var_temp=Mobj.h;
var_temp(test_i)=0;

%% Spheredist:

% Use spheredist and find nodes falling on the ib lines 
% Set h at those nodes = ice-bridge bathymetry

% There is no information available for the eastern half of the fjord. 
% We densify the IB lines (densify_iblines.m) assuming that the sill runs
% across the width of the fjord.

densify_iblines

ib_index=find(inpolygon(extrap_ib_lon,extrap_ib_lat,xp,yp)); %ice bridge lines inside N-S polygon
ib_bath_indexed=extrap_ib_bth(ib_index);
ib_lat_indexed=extrap_ib_lat(ib_index);ib_lon_indexed=extrap_ib_lon(ib_index);

for i=1:length(ib_lat_indexed)
    l=spheredist(Mobj.lat,Mobj.lon,ib_lat_indexed(i,1),ib_lon_indexed(i,1));
    ext_ind(i,1)=find(l==min(l));
    var_temp(ext_ind(i,1))=ib_bath_indexed(i,1);
end

%Visualise
figure(5);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 1200]);hold on;

%% Densify further:

% Note: var_temp was initialised with 0. Then, we populated it with IB
% bathymetry lines running across the width of the fjord. 

% Now, we want to add the Mobj.h values on the flanks.

% Make two polygons (one on each flank) and set h at the nodes inside 
% polygons to 'h' from Mobj (i.e. same as what we have from before)

% First, l = left

%[g_xp_l,g_yp_l] = ginput; 
%save g_xp_l; save g_yp_l;
load ./g_xp_l; load ./g_yp_l; 

%Visualise
figure(6);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 1200]);hold on;
hold on;
scatter(g_xp_l,g_yp_l,10,'k')

ind_g_p_l=find(inpolygon(Mobj.lon,Mobj.lat,g_xp_l,g_yp_l));
var_temp(ind_g_p_l)=Mobj.h(ind_g_p_l);
var_temp_lon_l=Mobj.lon(ind_g_p_l);var_temp_lat_l=Mobj.lat(ind_g_p_l);

%Visualise again
figure(7);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 1200]);hold on;
hold on;
scatter(g_xp_l,g_yp_l,10,'k')

% Then, r = right:
%[g_xp_r,g_yp_r] = ginput; 
%save g_xp_r; save g_yp_r;
load g_xp_r; load g_yp_r; 

%Visualise
figure(8);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 1200]);hold on;
hold on;
scatter(g_xp_r,g_yp_r,10,'k')

ind_g_p_r=find(inpolygon(Mobj.lon,Mobj.lat,g_xp_r,g_yp_r));
var_temp(ind_g_p_r)=Mobj.h(ind_g_p_r);
var_temp_lon_r=Mobj.lon(ind_g_p_r);var_temp_lat_r=Mobj.lat(ind_g_p_r);

%Visualise again
figure(9);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 1200]);hold on;
hold on;
scatter(g_xp_r,g_yp_r,10,'k')

% Note: var_temp should have IB + Mobj.h values at the desired nodes inside 
% our N-S polygon given by (temp_i). However, there seems to be a bug that 
% removes the added IB lines from var_temp. We can ignore this for now, as 
% we still have these variables logged separately (i.e. IB lines and Mobj.h
% polygons for the right and left flank):

% Create two new bathymetry variables (one for each flank)
var_h_r=Mobj.h(ind_g_p_r); var_h_l=Mobj.h(ind_g_p_l);

%% Get the Mobj.lat+lon+h values at the North and South crack lines:

crackline_N=NaN(size(Ext_Lon_Crack_N));

for i=1:length(Ext_Lon_Crack_N)
    l=spheredist(Mobj.lat,Mobj.lon,Ext_Lat_Crack_N(i,1),Ext_Lon_Crack_N(i,1));
    crackline_N(i,1)=find(l==min(l));
end

crackline_N_lat=Mobj.lat(crackline_N);crackline_N_lon=Mobj.lon(crackline_N);
crackline_N_h=Mobj.h(crackline_N);


crackline_S=NaN(size(Ext_Lon_Crack_S));

for i=1:length(Ext_Lon_Crack_S)
    l=spheredist(Mobj.lat,Mobj.lon,Ext_Lat_Crack_S(i,1),Ext_Lon_Crack_S(i,1));
    crackline_S(i,1)=find(l==min(l));
end

crackline_S_lat=Mobj.lat(crackline_S);crackline_S_lon=Mobj.lon(crackline_S);
crackline_S_h=Mobj.h(crackline_S);

%% Append everything together:

% IB Lines -> ib_lon_indexed/ib_lat_indexed/ib_bath_indexed
% Polygons -> var_temp_lon_r/var_temp_lat_r/var_h_r
%             var_temp_lon_l/var_temp_lat_l/var_h_l 

% &
% Crack -> crackline_N_lon,crackline_N_lat,crackline_N_h
%          crackline_S_lon,crackline_S_lat,crackline_S_h


interp_lons=[ib_lon_indexed;var_temp_lon_r;var_temp_lon_l;crackline_N_lon;crackline_S_lon];
interp_lats=[ib_lat_indexed;var_temp_lat_r;var_temp_lat_l;crackline_N_lat;crackline_S_lat];
interp_h=[ib_bath_indexed;var_h_r;var_h_l;crackline_N_h;crackline_S_h];


%% Interpolate:

% Gridded Interpolant:

% vq=griddata(interp_lons,interp_lats,interp_h,Mobj.lon(test_i),Mobj.lat(test_i),'natural');
% 
% figure(10);clf
% scatter(Mobj.lon(test_i),Mobj.lat(test_i),100,vq);colorbar;colormap jet;
% caxis([0 1200]);

% Scattered Interpolant:

F=scatteredInterpolant(interp_lons,interp_lats,interp_h,'natural');
h_val=F(Mobj.lon(test_i),Mobj.lat(test_i));

%Visualise:
figure(10);clf
scatter(Mobj.lon(test_i),Mobj.lat(test_i),10,h_val);colorbar;colormap jet;caxis([0 900]);


%% Test the final variable "new_h" before saving

new_h=Mobj.h;
new_h(test_i)=h_val;

figure(11);clf
plot_field(Mobj,new_h,'coordinate','spherical');colormap jet;caxis([0 900]);hold on;
scatter(Ext_Lon_Crack_N,Ext_Lat_Crack_N,10,'k');
scatter(Ext_Lon_Crack_S,Ext_Lat_Crack_S,10,'k');

% WCT:
old_wct=Mobj.h-Mobj.zisf;
new_wct=new_h-Mobj.zisf;

figure(12);clf
plot_field(Mobj,new_wct,'coordinate','spherical');colormap jet;caxis([0 800]);hold on;
scatter(Ext_Lon_Crack_N,Ext_Lat_Crack_N,10,'k');
scatter(Ext_Lon_Crack_S,Ext_Lat_Crack_S,10,'k');

figure(13);clf
plot_field(Mobj,old_wct,'coordinate','spherical');colormap jet;caxis([0 800]);hold on;
scatter(Ext_Lon_Crack_N,Ext_Lat_Crack_N,10,'k');
scatter(Ext_Lon_Crack_S,Ext_Lat_Crack_S,10,'k');

figure(14);clf
plot_field(Mobj,new_wct-old_wct,'coordinate','spherical');colormap jet;hold on;
scatter(Ext_Lon_Crack_N,Ext_Lat_Crack_N,10,'k');
scatter(Ext_Lon_Crack_S,Ext_Lat_Crack_S,10,'k');

%% Store in Mobj.h

Mobj.h=new_h;

%Visualise:
figure(8);clf
plot_field(Mobj,Mobj.h,'coordinate','spherical');colormap jet;caxis([0 900]);hold on;

%% Create new input files 

%Set up model specific components:
Mobj = setup_metrics(Mobj);

casename=get_cstr;

% Write for FVCOM:
% switch to sphercial coordinates
Mobj.nativeCoords='spherical';
% dump mesh and connectivity
write_FVCOM_grid(Mobj,'./Input_Files/Natural_Sill/grd.dat');
% dump bathymetry
write_FVCOM_bath(Mobj,'./Input_Files/Natural_Sill/dep.dat');
% % dump open boundary node list
write_FVCOM_obc(Mobj,'./Input_Files/Natural_Sill/obc.dat');
% dump sponge layer file
write_FVCOM_sponge(Mobj,'./Input_Files/Natural_Sill/spg.dat');
% dump Coriolis file
write_FVCOM_cor(Mobj,'./Input_Files/Natural_Sill/cor.dat');

fileprefix='pf_o_1';
write_FVCOM_iceshelfdraft(Mobj, fileprefix);
 
save ./Input_Files/Natural_Sill/IBN_M.mat Mobj;

%% Plot to check:

clear all;close all;clc;
load ./Input_Files/Natural_Sill/IBN_M.mat;

figure (1);clf
plot_field(Mobj,Mobj.zisf);caxis([0 300]);shading flat;colormap jet;

figure (2); clf
plot_field(Mobj,Mobj.h);caxis([0 800]);shading flat;colormap jet;

figure (3); clf
plot_field(Mobj,Mobj.h-Mobj.zisf);caxis([0 800]);shading flat;colormap jet;


%% Other tests:

% 1. Tried the smooth function: 
% z=smooth(new_h(indices_in_polygon),'method');

% 2. Applied Gaussian filter to smooth the data
% -> Create filter
% sigma = 10; % pick sigma value for the gaussian
% gaussFilter = gausswin(6*sigma + 1)';
% gaussFilter = gaussFilter / sum(gaussFilter); % normalize
% 
% filteredY = conv(new_h(test_i), gaussFilter, 'same');
% 
% new_h(test_i)=filteredY;
% 
% figure(12);clf
% plot_field(Mobj,new_h,'coordinate','spherical');colormap jet;