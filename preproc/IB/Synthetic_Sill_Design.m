%% Gaussian (Bell Curve) synthetic sill:

%%
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;

load ./M.mat
load ./Lon_Crack;load ./Lat_Crack;

%% Use the crack lines to create a polygon:

Ext_Lat_Crack_S1=Lat_Crack-0.02;Ext_Lon_Crack_S1=Lon_Crack+0.04;
Ext_Lat_Crack_S2=Lat_Crack-0.12;Ext_Lon_Crack_S2=Lon_Crack+0.32;

% Open the derived 'wct' figure to get a sense of where you want the sill:
openfig('Fig_Script.fig');hold on;
scatter(Lon_Crack,Lat_Crack,10,'k');
scatter(Ext_Lon_Crack_S1,Ext_Lat_Crack_S1,10,'k');
scatter(Ext_Lon_Crack_S2,Ext_Lat_Crack_S2,10,'r');

x_fill_l=-61.5;y_fill_l=80.75;
x_fill_r=-60.2;y_fill_r=80.9;

xp=[Ext_Lon_Crack_S1;x_fill_l;flipud(Ext_Lon_Crack_S2);x_fill_r];
yp=[Ext_Lat_Crack_S1;y_fill_l;flipud(Ext_Lat_Crack_S2);y_fill_r];


%Test the Polygon:
test_i =find(inpolygon(Mobj.lon,Mobj.lat,xp,yp));
var_temp=Mobj.h;
var_temp(test_i)=5000; %assign an absurd value

%Visualise
figure(2);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 1200]);hold on;
scatter(xp,yp,10,'k');

%% Create the 'sill' polygon:

%[g_xp_sill,g_yp_sill] = ginput; 
%save g_xp_sill;save g_yp_sill;
load g_xp_sill; load g_yp_sill;

ind_g_sill=find(inpolygon(Mobj.lon,Mobj.lat,g_xp_sill,g_yp_sill));
clear var_temp
var_temp=Mobj.h;
var_temp(test_i)=0;
var_temp(ind_g_sill)=350; % Chose the sill depth. IB = c.a. 450

%Visualise
figure(3);clf
plot_field(Mobj,var_temp,'coordinate','spherical');colormap jet;caxis([0 1200]);

% If everything looks okay, get the x,y and h for the sill:
xsill=Mobj.lon(ind_g_sill);ysill=Mobj.lat(ind_g_sill);hsill=var_temp(ind_g_sill);

%% Get the Mobj.lat+lon+h values at the two crack lines:

crackline_S1=NaN(size(Ext_Lon_Crack_S1));

for i=1:length(Ext_Lon_Crack_S1)
    l=spheredist(Mobj.lat,Mobj.lon,Ext_Lat_Crack_S1(i,1),Ext_Lon_Crack_S1(i,1));
    crackline_S1(i,1)=find(l==min(l));
end

crackline_S1_lat=Mobj.lat(crackline_S1);crackline_S1_lon=Mobj.lon(crackline_S1);
crackline_S1_h=Mobj.h(crackline_S1);

crackline_S2=NaN(size(Ext_Lon_Crack_S2));

for i=1:length(Ext_Lon_Crack_S2)
    l=spheredist(Mobj.lat,Mobj.lon,Ext_Lat_Crack_S2(i,1),Ext_Lon_Crack_S2(i,1));
    crackline_S2(i,1)=find(l==min(l));
end

crackline_S2_lat=Mobj.lat(crackline_S2);crackline_S2_lon=Mobj.lon(crackline_S2);
crackline_S2_h=Mobj.h(crackline_S2);

%% Interpolate:

% We have inside our polygon (test_i) -> x,y and h of the sill. x,y and h
% at the two crack lines (North and South bounds):

interp_lons=[xsill;crackline_S1_lon;crackline_S2_lon];
interp_lats=[ysill;crackline_S1_lat;crackline_S2_lat];
interp_h=[hsill;crackline_S1_h;crackline_S2_h];

F=scatteredInterpolant(interp_lons,interp_lats,interp_h,'natural');
h_val=F(Mobj.lon(test_i),Mobj.lat(test_i));

%Visualise:
figure(4);clf
scatter(Mobj.lon(test_i),Mobj.lat(test_i),10,h_val);colorbar;colormap jet;caxis([0 900]);

%% Get the new depth and visualise:

new_h=Mobj.h;
new_h(test_i)=h_val; %From interpolation

figure(5);clf
plot_field(Mobj,new_h,'coordinate','spherical');colormap jet;caxis([0 900]);hold on;

% WCT:
old_wct=Mobj.h-Mobj.zisf;
new_wct=new_h-Mobj.zisf;

figure(6);clf
plot_field(Mobj,old_wct,'coordinate','spherical');colormap jet;caxis([0 800]);hold on;

figure(7);clf
plot_field(Mobj,new_wct,'coordinate','spherical');colormap jet;caxis([0 800]);hold on;


%% Visualise the cross - section

%%%%%%%%%%%%%%%%%%%%% A section along the center %%%%%%%%%%%%%%%%%%%%%%%%%%
x1=299;
y1=80.88;
x2=299.38;
y2=80.72;
xn = 450;
xx = linspace(x1,x2,xn);
yy = linspace(y1,y2,xn);

%clear is2
for i = 1:numel(xx)
    dum = spheredist(yy(i),xx(i)-360,Mobj.lat,Mobj.lon-360);
    is(i) = find(dum == min(dum));
end
is = unique(is,'stable');

h_along_2=new_h(is);h_alonglat_2=Mobj.lat(is);h_alonglon_2=Mobj.lon(is);
h_along2_old=Mobj.h(is);

% For cross-section
zcrs = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-new_h),Mobj.siglay),Mobj.zisf); %change to old_h : use Mobj.h instead of new_h
zcrs = zcrs(is,:);
ll = repmat(h_alonglon_2,[1 23]); %or alonglat

dval=ones(size(ll)); %dummy fill values for pcolorjw

figure(100);clf
pcolorjw(ll,-zcrs,dval);
set(gca,'Color','k');
ylim([-1100 0]);

figure(101);clf
plot_field(Mobj,new_wct,'coordinate','spherical');colormap jet;caxis([0 800]);hold on;
scatter(x1-360,y1,140,'k','o','filled');
scatter(x2-360,y2,140,'k','o','filled');

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
write_FVCOM_grid(Mobj,'./Input_Files/Synthetic_Sill/IB_Exaggerated/grd.dat');
% dump bathymetry
write_FVCOM_bath(Mobj,'./Input_Files/Synthetic_Sill/IB_Exaggerated/dep.dat');
% % dump open boundary node list
write_FVCOM_obc(Mobj,'./Input_Files/Synthetic_Sill/IB_Exaggerated/obc.dat');
% dump sponge layer file
write_FVCOM_sponge(Mobj,'./Input_Files/Synthetic_Sill/IB_Exaggerated/spg.dat');
% dump Coriolis file
write_FVCOM_cor(Mobj,'./Input_Files/Synthetic_Sill/IB_Exaggerated/cor.dat');

fileprefix='pf_o_1';
write_FVCOM_iceshelfdraft(Mobj, fileprefix);

save ./Input_Files/Synthetic_Sill/IB_Exaggerated/IBE_M.mat Mobj;

%% Plot to check:

clear all;close all;clc;
load ./Input_Files/Synthetic_Sill/IB_Exaggerated/IBE_M.mat;

figure (1);clf
plot_field(Mobj,Mobj.zisf);caxis([0 300]);shading flat;colormap jet;

figure (2); clf
plot_field(Mobj,Mobj.h);caxis([0 800]);shading flat;colormap jet;

figure (3); clf
plot_field(Mobj,Mobj.h-Mobj.zisf);caxis([0 800]);shading flat;colormap jet;