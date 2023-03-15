%% Plotting Script:

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));

load ./M.mat; 

% Load all Atm-Forc variables:
load ./T_2m_air_forc ; load ./SW_forc; load ./LW_forc; 
load ./EV_forc; load ./RA_forc; load RH_forc; load P_Surf_forc;
load U10m_forc; load V10m_forc;

% Load ERA-RACMO boundary:
load ./rac_x_box.mat; load ./rac_y_box.mat;

%% Plot forcing outputs from a particular day of the model output:

%%%%% 1. Check the model output time_steps: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn='run1_uv_hf_fwf_0001.nc';
o_time=ncread(fn,'time'); o_time=o_time+678942;
o_dates=datestr(o_time);

%%%%% 2. Get the atmospheric forcing time-step(s) for that day: %%%%%%%%%%%

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/TA_RAC_2016.mat
atm_time=time; clear TA; clear time;

atm_dates=datestr(atm_time); % compare with o_dates

% The model out t_step = 10, where we are writing output once every 5 days, 
% corresponds to July 16 00:00:00 - or Day 46. In the atm forcing, this 
% corresponds to the atm t_step = 8 * 45 + 1 (8 3-hour t_steps each day)

plot_time=atm_time(361); 

%%%%%% 3. Get the atm forcing variables at the desired time step ("plot_time")

plot_2mair=T_2m_air_forc(:,361);plot_sw=SW_forc(:,361);

figure(1);clf
plot_field(Mobj,plot_2mair,'coordinate','spherical');colormap jet;
scatter(rac_x_box,rac_y_box,1,'k');
xlim([-120 20]);ylim([74 87]);

figure(2);clf
plot_field(Mobj,plot_sw,'coordinate','spherical');colormap jet;
scatter(rac_x_box,rac_y_box,1,'k');
xlim([-120 20]);ylim([74 87]);

%% Plot Time - Series of two nodes across the ERA-RACMO boundary:

%% 1. Get the nodes:
load ./Atm_Forc_T_Vec.mat

% From the scatter-plot (Atm_Forc_Outputs.m script), get coordinates:
x1=-74.01;y1=76.42; % Node 1: ERA
x2=-73.71;y2=76.46;  % Node 2: RACMO

% Find the nodes corresponding to the ERA and RACMO points
rac_dist=spheredist(Mobj.lat,Mobj.lon,y1,x1);
rac_ind=find(rac_dist==min(rac_dist));

era_dist=spheredist(Mobj.lat,Mobj.lon,y2,x2);
era_ind=find(era_dist==min(era_dist));

%% 2. Crop the forcing file at these nodes
SW_rac=SW_forc(rac_ind,:)';SW_era=SW_forc(era_ind,:)';
LW_rac=LW_forc(rac_ind,:)';LW_era=LW_forc(era_ind,:)';
EV_rac=EV_forc(rac_ind,:)';EV_era=EV_forc(era_ind,:)';
RA_rac=RA_forc(rac_ind,:)';RA_era=RA_forc(era_ind,:)';
RH_rac=RH_forc(rac_ind,:)';RH_era=RH_forc(era_ind,:)';
Psurf_rac=P_Surf_forc(rac_ind,:)';Psurf_era=P_Surf_forc(era_ind,:)';
U_rac=U10m_forc(rac_ind,:)';U_era=U10m_forc(era_ind,:)';
V_rac=V10m_forc(rac_ind,:)';V_era=V10m_forc(era_ind,:)';

%% 3.A) Plot for SW:

figure(1);clf
%Check where we are
%scatter(Mobj.lon,Mobj.lat,10,SW_forc(:,5));colorbar;colormap jet;hold on;
plot_field(Mobj,SW_forc(:,7),'coordinate','spherical');colormap jet;
hold on;
scatter(Mobj.lon(rac_ind),Mobj.lat(rac_ind),100,'*','r');
scatter(Mobj.lon(era_ind),Mobj.lat(era_ind),100,'*','k');

figure(2);clf %First 7 days from June 1 00:00:00 to June 8 00:00:00
plot(RACMO_time(1:57),SW_era(1:57),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(1:57),SW_rac(1:57),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('SW in W/m^2');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(3);clf %Last 7 days from August 24 00:00:00 to August 31 00:00:00
plot(RACMO_time(end-56:end),SW_era(end-56:end),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(end-56:end),SW_rac(end-56:end),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('SW in W/m^2');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(4);clf
plot(RACMO_time,SW_era,'r','LineWidth',2);hold on;
plot(RACMO_time,SW_rac,'k','LineWidth',1);
xlabel('Time'),ylabel('SW in W/m^2');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

%% Plot Other variables

figure(5);clf
plot(RACMO_time(1:57),LW_era(1:57),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(1:57),LW_rac(1:57),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('LW Down in W/m^2');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node','location','best');

figure(6);clf %Last 7 days from August 24 00:00:00 to August 31 00:00:00
plot(RACMO_time(end-56:end),LW_era(end-56:end),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(end-56:end),LW_rac(end-56:end),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('LW Down in W/m^2');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(7);clf
plot(RACMO_time(1:57),EV_era(1:57),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(1:57),EV_rac(1:57),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Evap in m/s');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node','location','best');

figure(8);clf %Last 7 days from August 24 00:00:00 to August 31 00:00:00
plot(RACMO_time(end-56:end),EV_era(end-56:end),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(end-56:end),EV_rac(end-56:end),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Evap in m/s');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(9);clf
plot(RACMO_time(1:57),RA_era(1:57),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(1:57),RA_rac(1:57),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Precip in m/s');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(10);clf %Last 7 days from August 24 00:00:00 to August 31 00:00:00
plot(RACMO_time(end-56:end),RA_era(end-56:end),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(end-56:end),RA_rac(end-56:end),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Precip in m/s');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(11);clf
plot(RACMO_time(1:57),RH_era(1:57),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(1:57),RH_rac(1:57),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('RH');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(12);clf %Last 7 days from August 24 00:00:00 to August 31 00:00:00
plot(RACMO_time(end-56:end),RH_era(end-56:end),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(end-56:end),RH_rac(end-56:end),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('RH');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(13);clf
plot(RACMO_time(1:57),Psurf_era(1:57),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(1:57),Psurf_rac(1:57),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Psurf in Pa');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(14);clf %Last 7 days from August 24 00:00:00 to August 31 00:00:00
plot(RACMO_time(end-56:end),Psurf_era(end-56:end),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(end-56:end),Psurf_rac(end-56:end),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Psurf in Pa');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(15);clf
plot(RACMO_time(1:57),U_era(1:57),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(1:57),U_rac(1:57),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Uwind in m/s');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node','location','best');

figure(16);clf %Last 7 days from August 24 00:00:00 to August 31 00:00:00
plot(RACMO_time(end-56:end),U_era(end-56:end),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(end-56:end),U_rac(end-56:end),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Uwind in m/s');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node','location','best');

figure(17);clf
plot(RACMO_time(1:57),V_era(1:57),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(1:57),V_rac(1:57),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Vwind in m/s');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');

figure(18);clf %Last 7 days from August 24 00:00:00 to August 31 00:00:00
plot(RACMO_time(end-56:end),V_era(end-56:end),'-rs','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','red');hold on;
plot(RACMO_time(end-56:end),V_rac(end-56:end),'-ko','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black');hold on;
xlabel('Time'),ylabel('Vwind in m/s');
datetick('x','ddmmmyyyy');
legend('racmo-node','era-node');