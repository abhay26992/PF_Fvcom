%% unm intensification time-series

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
%addpath(genpath('E:/Studies/Model_Out/Gamma_0.1/fvcom_4/Nest+Ice/PF_stdyr_runs/Output_years'));

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Bathymetry_IceBridge/Input_Files/Natural_Sill/IBN_M.mat

%% Seasons:

iw=30:105; %winter
is=170:245; %summer

%% Prepare vtvars

load ./mtransport/vtvars_std.mat; vt_std=vtvars; clear vtvars;
load ./mtransport/vtvars_qsa.mat; vt_qsa=vtvars; clear vtvars;
load ./mtransport/vtvars_qsb.mat; vt_qsb=vtvars; clear vtvars;

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

%% Transport calculation

% Theta: Angle between the successive section points (cell-indices) & x-axis

xc=Mobj.xc(vt_std.ind);yc=Mobj.yc(vt_std.ind);

for i=1:length(vt_std.ind)-1
    
    theta=atan2(yc(i+1)-yc(i),xc(i+1)-xc(i)); % breaking it into (#cell indices - 1 = 179) small sections
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Standard Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    unm_l=squeeze(vt_std.ux(i,:,:)).*sin(theta)-squeeze(vt_std.vy(i,:,:)).*cos(theta);
    unm_r=squeeze(vt_std.ux(i+1,:,:)).*sin(theta)-squeeze(vt_std.vy(i+1,:,:)).*cos(theta);
    
    
    for k=1:size(vt_std.ux,2)
        for l=1:size(vt_std.ux,3)
            
            if unm_l(k,l)>0
            unm_l(k,l)=0;
            end
            
            if unm_r(k,l)>0
            unm_r(k,l)=0;
            end
            
        end   
    end
    
    unm_std_av(i,:,:)=(unm_l+unm_r)/2; %avg. over successive small sections
    
    
    
    clear unm_l unm_r
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% qsa run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    unm_l=squeeze(vt_qsa.ux(i,:,:)).*sin(theta)-squeeze(vt_qsa.vy(i,:,:)).*cos(theta);
    unm_r=squeeze(vt_qsa.ux(i+1,:,:)).*sin(theta)-squeeze(vt_qsa.vy(i+1,:,:)).*cos(theta);
    
    
    for k=1:size(vt_qsa.ux,2)
        for l=1:size(vt_qsa.ux,3)
            
            if unm_l(k,l)>0
            unm_l(k,l)=0;
            end
            
            if unm_r(k,l)>0
            unm_r(k,l)=0;
            end
            
        end   
    end
    
    unm_qsa_av(i,:,:)=(unm_l+unm_r)/2; %avg. over successive small sections
        
    
    
   clear unm_l unm_r
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% qsb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   unm_l=squeeze(vt_qsb.ux(i,:,:)).*sin(theta)-squeeze(vt_qsb.vy(i,:,:)).*cos(theta);
   unm_r=squeeze(vt_qsb.ux(i+1,:,:)).*sin(theta)-squeeze(vt_qsb.vy(i+1,:,:)).*cos(theta);
    
    
    for k=1:size(vt_qsb.ux,2)
        for l=1:size(vt_qsb.ux,3)
            
            if unm_l(k,l)>0
            unm_l(k,l)=0;
            end
            
            if unm_r(k,l)>0
            unm_r(k,l)=0;
            end
            
        end   
    end
    
   unm_qsb_av(i,:,:)=(unm_l+unm_r)/2; %avg. over successive small sections
    
    
   clear unm_l unm_r
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %For plotting
    px(i,1)=(xc(i+1)+xc(i))/2; 
    py(i,1)=(yc(i+1)+yc(i))/2;
        
                
end

%% Prepare plot vars 

%1. Coordinate components
[plat,plon]=polarstereo_inv(repmat(px,[1 size(Mobj.siglay,2)]),repmat(py,[1 size(Mobj.siglay,2)]));
siglayc_ind=siglayc(vt_std.ind,:);
hc_ind=hc(vt_std.ind,:);%hc_ind=repmat(hc_ind, [1 size(Mobj.siglay,2)]);

zi=ncread(fn,'zisf');
zi_c=node2cell(zi,Mobj.tri);
zi_c_ind=zi_c(vt_std.ind,:);
zisf_c_ind = bsxfun(@plus,bsxfun(@times,(zi_c_ind-hc_ind),siglayc_ind),zi_c_ind);

for i=1:length(vt_std.ind)-1
    p_zisf_c_ind(i,:)=(zisf_c_ind(i+1,:)+zisf_c_ind(i,:))/2;
end

% 2. Variables

% 2a. Mean over our periods of interest (summer vs winter) of the velocity

unm_std_av_io_summer=mean(unm_std_av(:,:,is),3); 
unm_std_av_io_winter=mean(unm_std_av(:,:,iw),3); 

unm_qsa_av_io_summer=mean(unm_qsa_av(:,:,is),3); 
unm_qsa_av_io_winter=mean(unm_qsa_av(:,:,iw),3); 

unm_qsb_av_io_summer=mean(unm_qsb_av(:,:,is),3); 
unm_qsb_av_io_winter=mean(unm_qsb_av(:,:,iw),3);

%% Plot

%vfilt description:

%-> Using vfilt with zero-padding ('nonans'), and
%-> Filter size = 6 (total ~180 points over ~15 km. So, 12 points/km is the
%approximate lateral resolution. So, 6 points ~ 0.5 km)
%-> Not using a median filter - Stripey results seen even with zero padding

%% 1. Summer Anomalies:


figure(1);clf
yfilt_1=vfilt(unm_qsa_av_io_summer-unm_std_av_io_summer,8,'nonans');
pcolorjw(plon+360,-p_zisf_c_ind,yfilt_1);c=colorbar;set(gca,'Color','k');
xlabel('Longitude');ylabel('Depth [m]');c.Label.String='Normal Velocity [m/s]';c.Label.FontSize = 22;
caxis([-0.06 0.01]);
cmocean('balance','pivot',0);
% colormap(brewermap(64,'*BrBG'));
% namecpt='temp_diff_18lev.cpt';
% cptcmap(namecpt,'mapping','direct','ncol',10);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',22,'FontWeight','Bold', 'LineWidth', 2);
title('Outflow - Summer')

figure(2);clf
yfilt_2=vfilt(unm_qsb_av_io_summer-unm_qsa_av_io_summer,8,'nonans');
pcolorjw(plon+360,-p_zisf_c_ind,yfilt_2);c=colorbar;set(gca,'Color','k');
xlabel('Longitude');ylabel('Depth [m]');c.Label.String='Normal Velocity [m/s]';c.Label.FontSize = 22;
caxis([-0.015 0.005]);
cmocean('balance','pivot',0);
% colormap(brewermap(64,'*RdYlBu'));
% namecpt='temp_diff_18lev.cpt';
% cptcmap(namecpt,'mapping','direct','ncol',10);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',22,'FontWeight','Bold', 'LineWidth', 2);


%% 2. Winter Anomalies


figure(3);clf
yfilt_3=vfilt(unm_qsa_av_io_winter-unm_std_av_io_winter,8,'nonans');
pcolorjw(plon+360,-p_zisf_c_ind,yfilt_3);c=colorbar;set(gca,'Color','k');
xlabel('Longitude');ylabel('Depth [m]');c.Label.String='Normal Velocity [m/s]';c.Label.FontSize = 22;
caxis([-0.06 0.01]);
cmocean('balance','pivot',0);
% colormap(brewermap(64,'*RdYlBu'));
% namecpt='temp_diff_18lev.cpt';
% cptcmap(namecpt,'mapping','direct','ncol',4);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',22,'FontWeight','Bold', 'LineWidth', 2);
title('Outflow - Winter')

figure(4);clf
yfilt_4=vfilt(unm_qsb_av_io_winter-unm_qsa_av_io_winter,8,'nonans');
pcolorjw(plon+360,-p_zisf_c_ind,yfilt_4);c=colorbar;set(gca,'Color','k');
xlabel('Longitude');ylabel('Depth [m]');c.Label.String='Normal Velocity [m/s]';c.Label.FontSize = 22;
caxis([-0.015 0.005]);
cmocean('balance','pivot',0);
% colormap(brewermap(64,'*RdYlBu'));
% namecpt='temp_diff_18lev.cpt';
% cptcmap(namecpt,'mapping','direct','ncol',10);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',22,'FontWeight','Bold', 'LineWidth', 2);

