% Save the merged (RACMO and ERA) atmospheric forcing outputs

clear all;close all;clc;

load ./M.mat
fname='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/evap.nc';
lon=ncread(fname,'lon');lat=ncread(fname,'lat');

%% 2-m air temperature

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/TA_RAC_2016.mat
RACMO_out=TA; RACMO_time=time; clear TA; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/TA_ERA_2016.mat
ERA_out=TA_ERA; ERA_time=time_era;clear TA_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"node");


T_2m_air_forc=merged_out; clear merged_out;
save T_2m_air_forc T_2m_air_forc

figure(1);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,T_2m_air_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Surface Pressure
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/PA_RAC_2016.mat
RACMO_out=PA; RACMO_time=time; clear PA; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/PA_ERA_2016.mat
ERA_out=PA_ERA; ERA_time=time_era;clear PA_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"node");

P_Surf_forc=merged_out; clear merged_out;
save P_Surf_forc P_Surf_forc

figure(2);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,P_Surf_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Relative Humidity
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/RH_RAC_2016.mat
RACMO_out=RH; RACMO_time=time; clear RH; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/RH_ERA_2016.mat
ERA_out=RH_ERA; ERA_time=time_era;clear RH_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"node");

RH_forc=merged_out; clear merged_out;
save RH_forc RH_forc

figure(3);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,RH_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% U10m 
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/UW_RAC_2016.mat
RACMO_out=UW; RACMO_time=time; clear UW; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/UW_ERA_2016.mat
ERA_out=UW_ERA; ERA_time=time_era;clear UW_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"cell");

U10m_forc=merged_out; clear merged_out;
save U10m_forc U10m_forc

figure(4);clf
load rac_x_box 
load rac_y_box
[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
scatter(lonc,latc,10,U10m_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% V10m

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/VW_RAC_2016.mat
RACMO_out=VW; RACMO_time=time; clear VW; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/VW_ERA_2016.mat
ERA_out=VW_ERA; ERA_time=time_era;clear VW_ERA; clear time_era;

hs_ERA_out=zeros(length(ERA_out),length(RACMO_time));

for p=1:length(ERA_out)
    hs_ERA_out(p,:)=interp1(ERA_time',ERA_out(p,:),RACMO_time');
end

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,hs_ERA_out,"cell");

V10m_forc=merged_out; clear merged_out;
save V10m_forc V10m_forc

figure(5);clf
load rac_x_box 
load rac_y_box
[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
scatter(lonc,latc,10,V10m_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% LW Down

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/LW_RAC_2016.mat
RACMO_out=d_LW; RACMO_time=time; clear d_LW; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/LW_ERA_2016.mat
ERA_out=LW_ERA; ERA_time=time_era(2:end,1);clear LW_ERA; clear time_era;

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,"node");

LW_forc=merged_out; clear merged_out;
save LW_forc LW_forc

figure(6);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,LW_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% SW Net

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/SW_RAC_2016.mat
RACMO_out=SW; RACMO_time=time; clear SW; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/SW_ERA_2016.mat
ERA_out=SW_ERA; ERA_time=time_era;clear SW_ERA; clear time_era;

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,"node");

SW_forc=merged_out; clear merged_out;
save SW_forc SW_forc

figure(7);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,SW_forc(:,7));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Precipitation:

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/RA_RAC_2016.mat
RACMO_out=RA; RACMO_time=time; clear RA; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/TP_ERA_2016.mat
ERA_out=TP_ERA; ERA_time=time_era;clear TP_ERA; clear time_era;

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,"node");

RA_forc=merged_out; clear merged_out;
save RA_forc RA_forc

figure(8);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,RA_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Evaporation:

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/Evap_RAC_2016.mat
RACMO_out=Evap; RACMO_time=time; clear Evap; clear time;
load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/Evap_ERA_2016.mat
ERA_out=Evap_ERA; ERA_time=time_era;clear Evap_ERA; clear time_era;

merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,"node");

EV_forc=merged_out; clear merged_out;
save EV_forc EV_forc

figure(9);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,EV_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);

%% Merged P-E Map:

load ./RA_forc
load ./EV_forc

PE_forc= RA_forc - EV_forc;
save PE_forc PE_forc

figure(10);clf
load rac_x_box 
load rac_y_box
scatter(Mobj.lon,Mobj.lat,10,PE_forc(:,1));hold on;
scatter(rac_x_box,rac_y_box,1,'k');
shading flat;colorbar;colormap jet;
xlim([-120 20]);ylim([74 87]);