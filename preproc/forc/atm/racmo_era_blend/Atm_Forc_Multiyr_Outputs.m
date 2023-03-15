% Save the merged (RACMO and ERA) atmospheric forcing outputs

clear all;close all;clc;

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));

load ./M.mat
fname='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/evap.nc';
lon=ncread(fname,'lon');lat=ncread(fname,'lat');

%% 2-m air temperature

load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/TA_RAC_2016.mat
RACMO_out=TA; RACMO_time=time; clear TA; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/TA_ERA_2016.mat
ERA_out=TA_ERA; ERA_time=time_era;clear TA_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"node");


T_2m_air_forc=merged_out; clear merged_out;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/T_2m_air_forc_16 T_2m_air_forc -v7.3

figure(1);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,T_2m_air_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Surface Pressure
load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/PA_RAC_2016.mat
RACMO_out=PA; RACMO_time=time; clear PA; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/PA_ERA_2016.mat
ERA_out=PA_ERA; ERA_time=time_era;clear PA_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"node");

P_Surf_forc=merged_out; clear merged_out;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/P_Surf_forc_16 P_Surf_forc -v7.3

figure(2);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,P_Surf_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Relative Humidity
load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/RH_RAC_2016.mat
RACMO_out=RH; RACMO_time=time; clear RH; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/RH_ERA_2016.mat
ERA_out=RH_ERA; ERA_time=time_era;clear RH_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"node");

RH_forc=merged_out; clear merged_out;

save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/RH_forc_16 RH_forc -v7.3

figure(3);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,RH_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% U10m 
load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/UW_RAC_2016.mat
RACMO_out=UW; RACMO_time=time; clear UW; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/UW_ERA_2016.mat
ERA_out=UW_ERA; ERA_time=time_era;clear UW_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"cell");

U10m_forc=merged_out; clear merged_out;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/U10m_forc_16 U10m_forc -v7.3

figure(4);clf
load rac_x_box 
load rac_y_box
[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
scatter(lonc,latc,10,U10m_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% V10m

load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/VW_RAC_2016.mat
RACMO_out=VW; RACMO_time=time; clear VW; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/VW_ERA_2016.mat
ERA_out=VW_ERA; ERA_time=time_era;clear VW_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"cell");

V10m_forc=merged_out; clear merged_out;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/V10m_forc_16 V10m_forc -v7.3

figure(5);clf
load rac_x_box 
load rac_y_box
[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
scatter(lonc,latc,10,V10m_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% LW Down

load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/LW_RAC_2016.mat
RACMO_out=d_LW; RACMO_time=time; clear d_LW; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/LW_ERA_2016.mat
ERA_out=LW_ERA; ERA_time=time_era(2:end,1);clear LW_ERA; clear time_era;

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,"node");

LW_forc=merged_out; clear merged_out;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/LW_forc_16 LW_forc -v7.3

figure(6);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,LW_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% SW Net

load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/SW_RAC_2016.mat
RACMO_out=SW; RACMO_time=time; clear SW; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/SW_ERA_2016.mat
ERA_out=SW_ERA; ERA_time=time_era;clear SW_ERA; clear time_era;

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,"node");

SW_forc=merged_out; clear merged_out;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/SW_forc_16 SW_forc -v7.3

figure(7);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,SW_forc(:,7));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Precipitation:

load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/RA_RAC_2016.mat
RACMO_out=RA; RACMO_time=time; clear RA; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/TP_ERA_2016.mat
ERA_out=TP_ERA; ERA_time=time_era(2:end,1);clear TP_ERA; clear time_era;

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,"node");

RA_forc=merged_out; clear merged_out;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/RA_forc_16 RA_forc -v7.3

figure(8);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,RA_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Evaporation:

load D:/Studies/Atmospheric_Forcing/RACMO/RACMO_Forcing/2016/Evap_RAC_2016.mat
RACMO_out=Evap; RACMO_time=time; clear Evap; clear time;
load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/Evap_ERA_2016.mat
ERA_out=Evap_ERA; ERA_time=time_era(2:end,1);clear Evap_ERA; clear time_era;

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,"node");

EV_forc=merged_out; clear merged_out;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/EV_forc_16 EV_forc -v7.3

figure(9);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,EV_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Merged P-E Map:

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/RA_forc_16
load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/EV_forc_16

PE_forc= RA_forc - EV_forc;
save D:/Studies/Atmospheric_Forcing/Merged_Product/2016/PE_forc_16 PE_forc -v7.3

figure(10);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,PE_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);