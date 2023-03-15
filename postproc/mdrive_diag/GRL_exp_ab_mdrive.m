%% Melt rate and melt rate drivers

% Can be extended to include additional experiments

% Uses jlab

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));

%% Data

load('./mdrive/melt_std.mat')
load('./mdrive/melt_qsa.mat')
load('./mdrive/melt_qsb.mat')

load('./mdrive/Tdrive_std.mat')
load('./mdrive/Tdrive_qsa.mat')
load('./mdrive/Tdrive_qsb.mat')

load('./mdrive/ustar_std.mat')
load('./mdrive/ustar_qsa.mat')
load('./mdrive/ustar_qsb.mat')


load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Bathymetry_IceBridge/Input_Files/Natural_Sill/IBN_M.mat

zz=find(Mobj.zisf==0);
ustar_std(zz,:)=NaN;ustar_qsa(zz,:)=NaN;ustar_qsb(zz,:)=NaN;
Tdrive_std(zz,:)=NaN;Tdrive_qsa(zz,:)=NaN;Tdrive_qsb(zz,:)=NaN;
melt_std(zz,:)=NaN;melt_qsa(zz,:)=NaN;melt_qsb(zz,:)=NaN;


%% Mobj subset

regx =[-2.9 -2.5]*10^5;
regy =[-10 -9.1]*10^5;

nMobj = subset_mobj(Mobj,regx,regy);
nMobj.zisf=Mobj.zisf(nMobj.in_sub);
nMobj.lon=nMobj.lon+360; %To be consistent with other plots in the paper


%% Trim

cc=3; % <--- change if additional exps. are added.

% winter and summer months
iw = 30:105;
is = 170:245;

% TKL annual mean melt and seasonal anomalies wrt the annual mean:
secpera = 3600*24*365;

clear mm mas maw tdm tdas tdaw usdm

for r = 1:cc
    
    switch r             
       
            
        case 1
            
            run = 'std';
            mm(:,r) = vmean(melt_std(nMobj.in_sub,:),2)*secpera; % annual mean
            mas(:,r) = vmean(melt_std(nMobj.in_sub,is),2)*secpera-mm(:,r); % summer anomaly
            ms(:,r) = vmean(melt_std(nMobj.in_sub,is),2)*secpera;
            maw(:,r) = vmean(melt_std(nMobj.in_sub,iw),2)*secpera-mm(:,r); % win anomaly
            mw(:,r) = vmean(melt_std(nMobj.in_sub,iw),2)*secpera;

            % seasonal anomalies in thermal driving and friction velocity
            tdm(:,r) = vmean(Tdrive_std(nMobj.in_sub,:),2);
            tdas(:,r) = vmean(Tdrive_std(nMobj.in_sub,is),2)-tdm(:,r);
            tds(:,r) = vmean(Tdrive_std(nMobj.in_sub,is),2);
            tdaw(:,r) = vmean(Tdrive_std(nMobj.in_sub,iw),2)-tdm(:,r);
            tdw(:,r) = vmean(Tdrive_std(nMobj.in_sub,iw),2);

            usm(:,r) = vmean(ustar_std(nMobj.in_sub,:),2);
            usas(:,r) = vmean(ustar_std(nMobj.in_sub,is),2)-usm(:,r);
            uss(:,r) = vmean(ustar_std(nMobj.in_sub,is),2);
            usaw(:,r) = vmean(ustar_std(nMobj.in_sub,iw),2)-usm(:,r);
            usw(:,r) = vmean(ustar_std(nMobj.in_sub,iw),2);
            
        
            
        case 2
            
            run = 'qsa';
            mm(:,r) = vmean(melt_qsa(nMobj.in_sub,:),2)*secpera; % annual mean
            mas(:,r) = vmean(melt_qsa(nMobj.in_sub,is),2)*secpera-mm(:,r); % summer anomaly
            ms(:,r) = vmean(melt_qsa(nMobj.in_sub,is),2)*secpera;
            maw(:,r) = vmean(melt_qsa(nMobj.in_sub,iw),2)*secpera-mm(:,r); % win anomaly
            mw(:,r) = vmean(melt_qsa(nMobj.in_sub,iw),2)*secpera;

            % seasonal anomalies in thermal driving and friction velocity
            tdm(:,r) = vmean(Tdrive_qsa(nMobj.in_sub,:),2);
            tdas(:,r) = vmean(Tdrive_qsa(nMobj.in_sub,is),2)-tdm(:,r);
            tds(:,r) = vmean(Tdrive_qsa(nMobj.in_sub,is),2);
            tdaw(:,r) = vmean(Tdrive_qsa(nMobj.in_sub,iw),2)-tdm(:,r);
            tdw(:,r) = vmean(Tdrive_qsa(nMobj.in_sub,iw),2);

            usm(:,r) = vmean(ustar_qsa(nMobj.in_sub,:),2);
            usas(:,r) = vmean(ustar_qsa(nMobj.in_sub,is),2)-usm(:,r);
            uss(:,r) = vmean(ustar_qsa(nMobj.in_sub,is),2);
            usaw(:,r) = vmean(ustar_qsa(nMobj.in_sub,iw),2)-usm(:,r);
            usw(:,r) = vmean(ustar_qsa(nMobj.in_sub,iw),2);
            
        case 3
            
            run= 'qsb';
            mm(:,r) = vmean(melt_qsb(nMobj.in_sub,:),2)*secpera; % annual mean
            mas(:,r) = vmean(melt_qsb(nMobj.in_sub,is),2)*secpera-mm(:,r); % summer anomaly
            ms(:,r) = vmean(melt_qsb(nMobj.in_sub,is),2)*secpera;
            maw(:,r) = vmean(melt_qsb(nMobj.in_sub,iw),2)*secpera-mm(:,r); % win anomaly
            mw(:,r) = vmean(melt_qsb(nMobj.in_sub,iw),2)*secpera;

            % seasonal anomalies in thermal driving and friction velocity
            tdm(:,r) = vmean(Tdrive_qsb(nMobj.in_sub,:),2);
            tdas(:,r) = vmean(Tdrive_qsb(nMobj.in_sub,is),2)-tdm(:,r);
            tds(:,r) = vmean(Tdrive_qsb(nMobj.in_sub,is),2);
            tdaw(:,r) = vmean(Tdrive_qsb(nMobj.in_sub,iw),2)-tdm(:,r);
            tdw(:,r) = vmean(Tdrive_qsb(nMobj.in_sub,iw),2);

            usm(:,r) = vmean(ustar_qsb(nMobj.in_sub,:),2);
            usas(:,r) = vmean(ustar_qsb(nMobj.in_sub,is),2)-usm(:,r);
            uss(:,r) = vmean(ustar_qsb(nMobj.in_sub,is),2);
            usaw(:,r) = vmean(ustar_qsb(nMobj.in_sub,iw),2)-usm(:,r);
            usw(:,r) = vmean(ustar_qsb(nMobj.in_sub,iw),2);            
            
            
    end
    
    
end


%% Plots

%% 1a. Standard Run

figure(1);clf
plot_field(nMobj,mm(:,1),'coordinate','spherical');
caxis([0 100]);
colormap(brewermap(1000,'YlOrRd'));
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Longitude');ylabel('Latitude');

figure(2);clf
plot_field(nMobj,mas(:,1),'coordinate','spherical');
caxis([-30 30]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Longitude');ylabel('Latitude');

figure(3);clf
plot_field(nMobj,maw(:,1),'coordinate','spherical');
caxis([-30 30]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Longitude');ylabel('Latitude');

figure(4);clf
plot(nMobj.zisf,[tdaw(:,1),tdas(:,1)],'.');grid on;ylim([-0.2 0.2]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Ice Draft [m]');ylabel('\DeltaT [^{o}C]');
[objs, objh] = legend({'Winter','Summer'});
objhl = findobj(objh, 'type', 'line');
set(objhl, 'Markersize', 15);

figure(5);clf
plot_field(nMobj,tdas(:,1),'coordinate','spherical');
caxis(0.2*[-1 1]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Longitude');ylabel('Latitude');

figure(6);clf
plot_field(nMobj,tdaw(:,1),'coordinate','spherical');
caxis(0.2*[-1 1]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Longitude');ylabel('Latitude');

figure(7);clf
plot(nMobj.zisf,[usaw(:,1),usas(:,1)],'.');grid on;ylim([-4.5e-03 4.5e-03]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Ice Draft [m]');ylabel('u* [m/s]');
[objs, objh] = legend({'Winter','Summer'});
objhl = findobj(objh, 'type', 'line');
set(objhl, 'Markersize', 15);

figure(8);clf
plot_field(nMobj,usas(:,1),'coordinate','spherical');
caxis(0.004*[-1 1]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Longitude');ylabel('Latitude');

figure(9);clf
plot_field(nMobj,usaw(:,1),'coordinate','spherical');
caxis(0.004*[-1 1]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Longitude');ylabel('Latitude');
 

%% 1b. Qsa


figure(10);clf
subplot(3,3,1)
plot_field(nMobj,mm(:,2),'coordinate','spherical');
caxis([0 100]);
colormap(brewermap(1000,'YlOrRd'));
%c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(11);clf
subplot(3,3,2)
plot_field(nMobj,mas(:,2),'coordinate','spherical');
caxis([-30 30]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(12);clf
subplot(3,3,3)
plot_field(nMobj,maw(:,2),'coordinate','spherical');
caxis([-30 30]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(13);clf
subplot(3,3,4)
plot(nMobj.zisf,[tdaw(:,2),tdas(:,2)],'.');grid on;ylim([-0.2 0.2]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');
ylabel('\DeltaT [^{o}C]');
%[objs, objh] = legend({'Winter','Summer'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);

%figure(14);clf
subplot(3,3,5)
plot_field(nMobj,tdas(:,2),'coordinate','spherical');
caxis(0.2*[-1 1]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(15);clf
subplot(3,3,6)
plot_field(nMobj,tdaw(:,2),'coordinate','spherical');
caxis(0.2*[-1 1]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(16);clf
subplot(3,3,7)
plot(nMobj.zisf,[usaw(:,2),usas(:,2)],'.');grid on;ylim([-4.5e-03 4.5e-03]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Ice Draft [m]');ylabel('u* [m/s]');
%[objs, objh] = legend({'Winter','Summer'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);

%figure(17);clf
subplot(3,3,8)
plot_field(nMobj,usas(:,2),'coordinate','spherical');
caxis(0.004*[-1 1]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(18);clf
subplot(3,3,9)
plot_field(nMobj,usaw(:,2),'coordinate','spherical');
caxis(0.004*[-1 1]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');


%% 1c. Qsb


figure(19);clf
subplot(3,3,1)
plot_field(nMobj,mm(:,3),'coordinate','spherical');
caxis([0 100]);
colormap(brewermap(1000,'YlOrRd'));
%c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(20);clf
subplot(3,3,2)
plot_field(nMobj,mas(:,3),'coordinate','spherical');
caxis([-30 30]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(21);clf
subplot(3,3,3)
plot_field(nMobj,maw(:,3),'coordinate','spherical');
caxis([-30 30]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(22);clf
subplot(3,3,4)
plot(nMobj.zisf,[tdaw(:,3),tdas(:,3)],'.');grid on;ylim([-0.2 0.2]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');
ylabel('\DeltaT [^{o}C]');
%[objs, objh] = legend({'Winter','Summer'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);

%figure(23);clf
subplot(3,3,5)
plot_field(nMobj,tdas(:,3),'coordinate','spherical');
caxis(0.2*[-1 1]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(24);clf
subplot(3,3,6)
plot_field(nMobj,tdaw(:,3),'coordinate','spherical');
caxis(0.2*[-1 1]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(25);clf
subplot(3,3,7)
plot(nMobj.zisf,[usaw(:,3),usas(:,3)],'.');grid on;ylim([-4.5e-03 4.5e-03]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Ice Draft [m]');ylabel('u* [m/s]');
%[objs, objh] = legend({'Winter','Summer'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);

%figure(26);clf
subplot(3,3,8)
plot_field(nMobj,usas(:,3),'coordinate','spherical');
caxis(0.004*[-1 1]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(27);clf
subplot(3,3,9)
plot_field(nMobj,usaw(:,3),'coordinate','spherical');
caxis(0.004*[-1 1]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');


%% 2a. Melt Map anomalies for the different runs


figure(28);clf
subplot(3,3,1)
plot_field(nMobj,mm(:,1),'coordinate','spherical');
caxis([0 100]);
colormap(brewermap(1000,'YlOrRd'));
%c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(29);clf
subplot(3,3,2)
plot_field(nMobj,mm(:,2)-mm(:,1),'coordinate','spherical');
caxis([-25 125]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(30);clf
subplot(3,3,3)
plot_field(nMobj,mm(:,3)-mm(:,1),'coordinate','spherical');
caxis([-25 125]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(31);clf
subplot(3,3,4)
plot(nMobj.zisf,ms(:,2)-ms(:,1),'.');grid on;ylim([-100 500]);hold on;
plot(nMobj.zisf,ms(:,3)-ms(:,1),'.');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');ylabel('Melt rate [m/yr]');
%[objs, objh] = legend({'qsa-std','qsb-qsa'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);

%figure(32);clf
subplot(3,3,5)
plot_field(nMobj,ms(:,2)-ms(:,1),'coordinate','spherical');
caxis([-25 125]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(33);clf
subplot(3,3,6)
plot_field(nMobj,ms(:,3)-ms(:,1),'coordinate','spherical')
caxis([-25 125]);
cmocean('balance','pivot',0);
%c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(34);clf
subplot(3,3,7)
plot(nMobj.zisf,mw(:,2)-mw(:,1),'.');grid on;ylim([-100 500]);hold on;
plot(nMobj.zisf,mw(:,3)-mw(:,1),'.');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Ice Draft [m]');ylabel('Melt rate [m/yr]');
%[objs, objh] = legend({'qsa-std','qsb-qsa'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);

%figure(35);clf
subplot(3,3,8)
plot_field(nMobj,mw(:,2)-mw(:,1),'coordinate','spherical')
caxis([-25 125]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

%figure(36);clf
subplot(3,3,9)
plot_field(nMobj,mw(:,3)-mw(:,1),'coordinate','spherical')
caxis([-25 125]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'Melt rate [m/yr]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');


%% 2b. Plot_Field anomalies - Td


figure(37);clf
subplot(3,2,1)
plot_field(nMobj,tdm(:,2)-tdm(:,1),'coordinate','spherical');
caxis([-0.2 0.45]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

subplot(3,2,2)
plot_field(nMobj,tdm(:,3)-tdm(:,1),'coordinate','spherical');
%caxis([-0.05 0.05]); - this was the caxis for 3-2
caxis([-0.2 0.45]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

% figure(1000);clf
% plot_field(nMobj,tds(:,1),'coordinate','spherical');caxis([0.05 1.15]);
% cmocean('thermal');
% c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
% xlabel('Longitude');ylabel('Latitude');

subplot(3,2,3)
plot_field(nMobj,tds(:,2)-tds(:,1),'coordinate','spherical');
caxis([-0.2 0.45]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

subplot(3,2,4)
plot_field(nMobj,tds(:,3)-tds(:,1),'coordinate','spherical');
caxis([-0.2 0.45]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

subplot(3,2,5)
plot_field(nMobj,tdw(:,2)-tdw(:,1),'coordinate','spherical');
caxis([-0.2 0.45]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

subplot(3,2,6)
plot_field(nMobj,tdw(:,3)-tdw(:,1),'coordinate','spherical');
caxis([-0.2 0.45]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = '\DeltaT [^{o}C]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');


%% 2c. Plot_Field anomalies - Ustar


figure(38);clf
subplot(3,2,1)
plot_field(nMobj,usm(:,2)-usm(:,1),'coordinate','spherical');
caxis(10^-3*[-2.5 25]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

subplot(3,2,2)
plot_field(nMobj,usm(:,3)-usm(:,1),'coordinate','spherical');
caxis(10^-3*[-2.5 25]);
%caxis(10^-3*[-3.5 3.5]); - this was the caxis for 3-2
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

% figure(1000);clf
% plot_field(nMobj,uss(:,1),'coordinate','spherical');caxis(10^-3*[1 17]);
% cmocean('speed');
% c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
% xlabel('Longitude');ylabel('Latitude');

subplot(3,2,3)
plot_field(nMobj,uss(:,2)-uss(:,1),'coordinate','spherical');
caxis(10^-3*[-2.5 25]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

subplot(3,2,4)
plot_field(nMobj,uss(:,3)-uss(:,1),'coordinate','spherical');
caxis(10^-3*[-2.5 25]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

subplot(3,2,5)
plot_field(nMobj,usw(:,2)-usw(:,1),'coordinate','spherical');
caxis(10^-3*[-2.5 25]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');

subplot(3,2,6)
plot_field(nMobj,usw(:,3)-usw(:,1),'coordinate','spherical');
caxis(10^-3*[-2.5 25]);
cmocean('balance','pivot',0);
c=colorbar;c.Label.String = 'u* [m/s]';c.Label.FontSize = 28;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Longitude');ylabel('Latitude');


%% 2d. Scatter anomalies (Melt, Td, Ustar)


figure(39);clf
subplot(3,3,1)
plot(nMobj.zisf,[mm(:,2)-mm(:,1),mm(:,3)-mm(:,2)],'.');grid on;ylim([-50 200]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Ice Draft [m]');ylabel('Melt rate [m/yr]');
title('Annual');

subplot(3,3,2)
plot(nMobj.zisf,[ms(:,2)-ms(:,1),ms(:,3)-ms(:,2)],'.');grid on;ylim([-50 200]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');ylabel('Melt rate [m/yr]');
title('Summer');

subplot(3,3,3)
plot(nMobj.zisf,[mw(:,2)-mw(:,1),mw(:,3)-mw(:,2)],'.');grid on;ylim([-50 200]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');ylabel('Melt rate [m/yr]');
%[objs, objh] = legend({'qsa-std','qsb-qsa'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);
title('Winter');

subplot(3,3,4)
plot(nMobj.zisf,[tdm(:,2)-tdm(:,1),tdm(:,3)-tdm(:,2)],'.');grid on;ylim([-0.2 0.45]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Ice Draft [m]');ylabel('\DeltaT [^{o}C]');

subplot(3,3,5)
plot(nMobj.zisf,[tds(:,2)-tds(:,1),tds(:,3)-tds(:,2)],'.');grid on;ylim([-0.2 0.45]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');ylabel('\DeltaT [^{o}C]');

subplot(3,3,6)
plot(nMobj.zisf,[tdw(:,2)-tdw(:,1),tdw(:,3)-tdw(:,2)],'.');grid on;ylim([-0.2 0.45]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',14,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');ylabel('\DeltaT [^{o}C]');
%[objs, objh] = legend({'qsa-std','qsb-qsa'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);

subplot(3,3,7)
plot(nMobj.zisf,[usm(:,2)-usm(:,1),usm(:,3)-usm(:,2)],'.');grid on;ylim(10^-3*[-2.5 25]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
xlabel('Ice Draft [m]');ylabel('u* [m/s]');

subplot(3,3,8)
plot(nMobj.zisf,[uss(:,2)-uss(:,1),uss(:,3)-uss(:,2)],'.');grid on;ylim(10^-3*[-2.5 25]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');ylabel('u* [m/s]');

subplot(3,3,9)
plot(nMobj.zisf,[usw(:,2)-usw(:,1),usw(:,3)-usw(:,2)],'.');grid on;ylim(10^-3*[-2.5 25]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);
%xlabel('Ice Draft [m]');ylabel('u* [m/s]');
%[objs, objh] = legend({'qsa-std','qsb-qsa'});
%objhl = findobj(objh, 'type', 'line');
%set(objhl, 'Markersize', 15);

