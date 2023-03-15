%% Overturning vertical profile

% Plots the net flow through each depth level (+ve = inflow, -ve = outflow)

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
%addpath(genpath('E:/Studies/Model_Out/Gamma_0.1/fvcom_4/Nest+Ice/PF_stdyr_runs/Output_years'));

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Bathymetry_IceBridge/Input_Files/Natural_Sill/IBN_M.mat

%% Data

load ./mtransport/vtvars_std.mat; vt_std=vtvars; clear vtvars;
load ./mtransport/vtvars_qsa.mat; vt_qsa=vtvars; clear vtvars;
load ./mtransport/vtvars_qsb.mat; vt_qsb=vtvars; clear vtvars;
load ./mtransport/vtvars_qsc.mat; vt_qsc=vtvars; clear vtvars;
load ./mtransport/vtvars_qsd.mat; vt_qsd=vtvars; clear vtvars;

% Seasonal means
iw=30:105;
is=170:245;

%% Center Variables 

% Load the 'center' variables (siglay center, siglev center and h_center) 
% from a restart file:

cwdir=pwd;
cd C:\PhD\FVCOM\Matlab_Repository\Petermann_Bathy\Output\Gamma_0.1\fvcom_4\Nest+Atm+Ice\PF_stdyr_runs\Transport_Diagnostics;
fn='ocnice_2016_restart_0367.nc';
siglayc=ncread(fn,'siglay_center');
siglevc=ncread(fn,'siglev_center');
hc=ncread(fn,'h_center');
cd(cwdir);

%For plotting, to be modified later ...
siglayc_z=siglayc(vt_std.ind,:).*hc(vt_std.ind); 

% Get depth values at siglev center at the transect indices:
siglevc_z=siglevc(vt_std.ind,:).*hc(vt_std.ind); 

% Get the siglevc layer thickness 
for i=1:length(vt_std.ind)
    siggdiff(i,:)=abs(diff(siglevc_z(i,:))); 
end

siggdiff=repmat(siggdiff,[1 1 size(vt_std.u,3)]);

% Transport for each layer (for all 180 transect indices & all 731 t-steps)

ut_dz_std=vt_std.ux.*siggdiff; 
vt_dz_std=vt_std.vy.*siggdiff;

ut_dz_qsa=vt_qsa.ux.*siggdiff; 
vt_dz_qsa=vt_qsa.vy.*siggdiff;

ut_dz_qsb=vt_qsb.ux.*siggdiff; 
vt_dz_qsb=vt_qsb.vy.*siggdiff;

ut_dz_qsc=vt_qsc.ux.*siggdiff; 
vt_dz_qsc=vt_qsc.vy.*siggdiff;

ut_dz_qsd=vt_qsd.ux.*siggdiff; 
vt_dz_qsd=vt_qsd.vy.*siggdiff;

%% Transport Calculation

% Theta: Angle between the successive section points (cell-indices) & x-axis
% Distance: Distance between the successive small sections

xc=Mobj.xc(vt_std.ind);yc=Mobj.yc(vt_std.ind);

for i=1:length(vt_std.ind)-1
    
    theta=atan2(yc(i+1)-yc(i),xc(i+1)-xc(i)); % breaking it into (#cell indices - 1 = 179) small sections

    for j=1:size(siglayc,2)
        
        unm_l_std=ut_dz_std(i,j,:)*sin(theta)-vt_dz_std(i,j,:)*cos(theta);
        unm_r_std=ut_dz_std(i+1,j,:)*sin(theta)-vt_dz_std(i+1,j,:)*cos(theta);
        unm_av_std(i,j,:)=(unm_l_std+unm_r_std)/2; %avg. over successive small sections
        
        unm_l_qsa=ut_dz_qsa(i,j,:)*sin(theta)-vt_dz_qsa(i,j,:)*cos(theta);
        unm_r_qsa=ut_dz_qsa(i+1,j,:)*sin(theta)-vt_dz_qsa(i+1,j,:)*cos(theta);
        unm_av_qsa(i,j,:)=(unm_l_qsa+unm_r_qsa)/2; %avg. over successive small sections
                
        unm_l_qsb=ut_dz_qsb(i,j,:)*sin(theta)-vt_dz_qsb(i,j,:)*cos(theta);
        unm_r_qsb=ut_dz_qsb(i+1,j,:)*sin(theta)-vt_dz_qsb(i+1,j,:)*cos(theta);
        unm_av_qsb(i,j,:)=(unm_l_qsb+unm_r_qsb)/2; %avg. over successive small sections
        
        unm_l_qsc=ut_dz_qsc(i,j,:)*sin(theta)-vt_dz_qsc(i,j,:)*cos(theta);
        unm_r_qsc=ut_dz_qsc(i+1,j,:)*sin(theta)-vt_dz_qsc(i+1,j,:)*cos(theta);
        unm_av_qsc(i,j,:)=(unm_l_qsc+unm_r_qsc)/2; %avg. over successive small sections
        
        unm_l_qsd=ut_dz_qsd(i,j,:)*sin(theta)-vt_dz_qsd(i,j,:)*cos(theta);
        unm_r_qsd=ut_dz_qsd(i+1,j,:)*sin(theta)-vt_dz_qsd(i+1,j,:)*cos(theta);
        unm_av_qsd(i,j,:)=(unm_l_qsd+unm_r_qsd)/2; %avg. over successive small sections
        
        siglayc_z_unmav(i,j)=(siglayc_z(i,j)+siglayc_z(i+1,j))/2;
        
        dist=sqrt((xc(i+1)-xc(i)).^2+(yc(i+1)-yc(i)).^2);
        tp_dx_std(i,j,:)=unm_av_std(i,j,:)*dist;%transport through each small section
        tp_dx_qsa(i,j,:)=unm_av_qsa(i,j,:)*dist;%transport through each small section  
        tp_dx_qsb(i,j,:)=unm_av_qsb(i,j,:)*dist;%transport through each small section
        tp_dx_qsc(i,j,:)=unm_av_qsc(i,j,:)*dist;%transport through each small section
        tp_dx_qsd(i,j,:)=unm_av_qsd(i,j,:)*dist;%transport through each small section
        
    end
    
    
end

plot_z=mean(siglayc_z_unmav,1);

vol_tp_std=sum(tp_dx_std,1); %in m3/s
vol_tp_sv_std=squeeze(vol_tp_std)/1e06; %in Sv

vol_tp_qsa=sum(tp_dx_qsa,1); %in m3/s
vol_tp_sv_qsa=squeeze(vol_tp_qsa)/1e06; %in Sv

vol_tp_qsb=sum(tp_dx_qsb,1); %in m3/s
vol_tp_sv_qsb=squeeze(vol_tp_qsb)/1e06; %in Sv

vol_tp_qsc=sum(tp_dx_qsc,1); %in m3/s
vol_tp_sv_qsc=squeeze(vol_tp_qsc)/1e06; %in Sv

vol_tp_qsd=sum(tp_dx_qsd,1); %in m3/s
vol_tp_sv_qsd=squeeze(vol_tp_qsd)/1e06; %in Sv

%% Plots

%%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare plot vars %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ovt_std=mean(vol_tp_sv_std,2);
ovt_std_is=mean(vol_tp_sv_std(:,is),2);
ovt_std_iw=mean(vol_tp_sv_std(:,iw),2);

ovt_qsa=mean(vol_tp_sv_qsa,2);
ovt_qsa_is=mean(vol_tp_sv_qsa(:,is),2);
ovt_qsa_iw=mean(vol_tp_sv_qsa(:,iw),2);

ovt_qsb=mean(vol_tp_sv_qsb,2);
ovt_qsb_is=mean(vol_tp_sv_qsb(:,is),2);
ovt_qsb_iw=mean(vol_tp_sv_qsb(:,iw),2);

ovt_qsc=mean(vol_tp_sv_qsc,2);
ovt_qsc_is=mean(vol_tp_sv_qsc(:,is),2);
ovt_qsc_iw=mean(vol_tp_sv_qsc(:,iw),2);

ovt_qsd=mean(vol_tp_sv_qsd,2);
ovt_qsd_is=mean(vol_tp_sv_qsd(:,is),2);
ovt_qsd_iw=mean(vol_tp_sv_qsd(:,iw),2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1);clf
plot(ovt_qsa,plot_z,'m','LineWidth',1.2);hold on;
plot(ovt_qsa_iw,plot_z,'b','LineWidth',1.2);
plot(ovt_qsa_is,plot_z,'r','LineWidth',1.2);
grid on;xlabel('Overturning [Sv]');ylabel('Depth');
legend('qsa-annual','qsa-winter','qsa-summer');

figure(2);clf
plot(ovt_qsb,plot_z,'m','LineWidth',1.2);hold on;
plot(ovt_qsb_iw,plot_z,'b','LineWidth',1.2);
plot(ovt_qsb_is,plot_z,'r','LineWidth',1.2);
grid on;xlabel('Overturning [Sv]');ylabel('Depth');
legend('qsb-annual','qsb-winter','qsb-summer');

figure(3);clf
plot(ovt_std_is,plot_z,':r','LineWidth',1.2);hold on;
plot(ovt_std_iw,plot_z,':b','LineWidth',1.2);

plot(ovt_qsa_is,plot_z,'--r','LineWidth',1.2);
plot(ovt_qsa_iw,plot_z,'--b','LineWidth',1.2);

plot(ovt_qsb_is,plot_z,'r','LineWidth',1.2);
plot(ovt_qsa_iw,plot_z,'b','LineWidth',1.2);

grid on;xlabel('Overturning [Sv]');ylabel('Depth');
legend('std-summer','std-winter','qsa-summer','qsa-winter','qsb-summer','qsb-winter');


%% Summer - all runs in one plot for Kappa

figure(4);clf
plot(ovt_std_is,plot_z,'c','LineWidth',1.6);hold on;
plot(ovt_qsa_is,plot_z,'--b','LineWidth',1.6);
plot(ovt_qsb_is,plot_z,':b','LineWidth',1.6);
plot(ovt_qsc_is,plot_z,'m','LineWidth',1.6);
plot(ovt_qsd_is,plot_z,'r','LineWidth',1.6);

grid on;xlabel('Overturning [Sv]');ylabel('Depth [m]');
legend('std-summer','sp_a','sp_b','tp_a','tp_b');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
