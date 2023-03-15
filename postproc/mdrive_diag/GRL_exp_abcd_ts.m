%% TS-Section plots (Discharge experiments)

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
%addpath(genpath('E:/Studies/Model_Out/Gamma_0.1/fvcom_4/Nest+Ice/PF_stdyr_runs'));

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Bathymetry_IceBridge/Input_Files/Natural_Sill/IBN_M.mat

%% Output info

iw=30:105; 
isum=170:245; 

fn_std='D:/Abhay/Model_Out/PF_std_yr_2016/ocnice_2016_avg_0001.nc';
fn_qsa='D:/Abhay/Model_Out/PF_Q_sg.exp_a/exp_a_avg_0001.nc';
fn_qsb='D:/Abhay/Model_Out/PF_Q_sg.exp_b/exp_b_avg_0001.nc';

fn_qsc='D:/Abhay/Model_Out/PF_Q_sg.exp_c/exp_c_avg_0001.nc';
fn_qsd='D:/Abhay/Model_Out/PF_Q_sg.exp_d/exp_d_avg_0001.nc';

%% Read vars:

si = ncread(fn_std,'siglay');

lon = ncread(fn_std,'lon');
lat = ncread(fn_std,'lat');

zi = ncread(fn_std,'zisf');
zi = zi(:,1);
h  = ncread(fn_std,'h');

t = ncread(fn_std,'temp');
s = ncread(fn_std,'salinity');

ta = ncread(fn_qsa,'temp');
sa = ncread(fn_qsa,'salinity');

tb = ncread(fn_qsb,'temp');
sb = ncread(fn_qsb,'salinity');

tc = ncread(fn_qsc,'temp');
sc = ncread(fn_qsc,'salinity');

td = ncread(fn_qsd,'temp');
sd = ncread(fn_qsd,'salinity');

%% Create a ROI:

%Fjord
xp = [295 302 302 295];
yp = [80 80 82 82];

ii_p = find(inpolygon(lon,lat,xp,yp));
z_p = bsxfun(@plus,bsxfun(@times,(zi-h),si),zi);

figure(1);clf
plot_field(Mobj,Mobj.h,'coordinate','spherical');colormap jet;caxis([0 900]);
hold on;scatter(Mobj.lon(ii_p),Mobj.lat(ii_p),10,'k');

%% Make Sections: 

% 1. In the fjord

x1_p=300;  %Pfjord-centre
y1_p=80.56;  %Pfjord-centre
x2_p=297.11; %Pfjord-centre
y2_p=81.4;   %Pfjord-centre

% x1_p = 360-60.18;  %Pfjord-west
% y1_p = 80.53;      %Pfjord-west
% x2_p = 360-63.86;  %Pfjord-west
% y2_p = 81.62;      %Pfjord-west

%% 3. Compute

xn=6000;
xx_p = linspace(x1_p,x2_p,xn);
yy_p = linspace(y1_p,y2_p,xn);

for i = 1:numel(xx_p)
    dum_p = spheredist(yy_p(i),xx_p(i)-360,lat(ii_p),lon(ii_p)-360);
    is_p(i) = find(dum_p == min(dum_p));
end
is_p = unique(is_p,'stable');

%plot(xx_p-360,yy_p,'r','LineWidth',2);hold off;


%% 4. Pfjord section plot:

%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. Winter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_winter_pf=mean(t(ii_p(is_p),:,iw),3);
s_winter_pf=mean(s(ii_p(is_p),:,iw),3);

t_winter_pfa=mean(ta(ii_p(is_p),:,iw),3);
s_winter_pfa=mean(sa(ii_p(is_p),:,iw),3);

t_winter_pfb=mean(tb(ii_p(is_p),:,iw),3);
s_winter_pfb=mean(sb(ii_p(is_p),:,iw),3);

t_winter_pfc=mean(tc(ii_p(is_p),:,iw),3);
s_winter_pfc=mean(sc(ii_p(is_p),:,iw),3);

t_winter_pfd=mean(td(ii_p(is_p),:,iw),3);
s_winter_pfd=mean(sd(ii_p(is_p),:,iw),3);

%%%%%%%%%%%%%%%%%%%%%%%%%%% 2. Summer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_summer_pf=mean(t(ii_p(is_p),:,isum),3);
s_summer_pf=mean(s(ii_p(is_p),:,isum),3);

t_summer_pfa=mean(ta(ii_p(is_p),:,isum),3);
s_summer_pfa=mean(sa(ii_p(is_p),:,isum),3);

t_summer_pfb=mean(tb(ii_p(is_p),:,isum),3);
s_summer_pfb=mean(sb(ii_p(is_p),:,isum),3);

t_summer_pfc=mean(tc(ii_p(is_p),:,isum),3);
s_summer_pfc=mean(sc(ii_p(is_p),:,isum),3);

t_summer_pfd=mean(td(ii_p(is_p),:,isum),3);
s_summer_pfd=mean(sd(ii_p(is_p),:,isum),3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3. Annual %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_annual_pf=mean(t(ii_p(is_p),:,:),3);
s_annual_pf=mean(s(ii_p(is_p),:,:),3);

t_annual_pfa=mean(ta(ii_p(is_p),:,:),3);
s_annual_pfa=mean(sa(ii_p(is_p),:,:),3);

t_annual_pfb=mean(tb(ii_p(is_p),:,:),3);
s_annual_pfb=mean(sb(ii_p(is_p),:,:),3);

t_annual_pfc=mean(tc(ii_p(is_p),:,:),3);
s_annual_pfc=mean(sc(ii_p(is_p),:,:),3);

t_annual_pfd=mean(td(ii_p(is_p),:,:),3);
s_annual_pfd=mean(sd(ii_p(is_p),:,:),3);

%% 5. Plots - PFjord section

%% 1st plot - for Kappa - Summer absolute temperature from std run

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
dens_c_summer=sw_dens0(s_summer_pf,t_summer_pf); %in kg/m^3 %choose exp here- exp a/b
figure(1);clf
cont_lvl_up=linspace(26.9,27.1,2);%set contour level
cont_lvl_down=linspace(27.1,27.5,5);%set contour level
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),t_summer_pf); %choose exp here- exp a/b
%hold on;plot(lon(ii_p(is_p)),-h(ii_p(is_p)));plot(lon(ii_p(is_p)),-zi(ii_p(is_p))); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;
caxis([-1.8 -0.4]);
c.Label.String = 'Temperature [{\circ}C]';c.Label.FontSize = 18;
cmocean('thermal');
xlabel('Longitude');ylabel('Depth [m]');
hold on;
[ccont,hcont]=contour(ll_custom_p,-z_p(ii_p(is_p),:),dens_c_summer-1000,cont_lvl_up,'ShowText','on','LineColor','k','LineWidth',0.2);
hcont.LevelList=round(hcont.LevelList,1);
clabel(ccont,hcont,'labelspacing',10.0e+100,'FontWeight','bold','Color','w','FontSize',26);
[ccont,hcont]=contour(ll_custom_p,-z_p(ii_p(is_p),:),dens_c_summer-1000,cont_lvl_down,'ShowText','on','LineColor','k','LineWidth',0.2);
hcont.LevelList=round(hcont.LevelList,1);
clabel(ccont,hcont,'labelspacing',10.0e+10,'FontWeight','bold','Color','w','FontSize',26);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',20,'FontWeight','Bold', 'LineWidth', 2);


%% 5a. Annual absolute - Qsg run - Experiment A/B


ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
dens_c_annual=sw_dens0(s_annual_pfa,t_annual_pfa); %in kg/m^3 %choose exp here- exp a/b
figure(1);clf
cont_lvl_up=linspace(26.9,27.1,2);%set contour level
cont_lvl_down=linspace(27.1,27.5,5);%set contour level
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),t_annual_pfa); %choose exp here- exp a/b
%hold on;plot(lon(ii_p(is_p)),-h(ii_p(is_p)));plot(lon(ii_p(is_p)),-zi(ii_p(is_p))); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;
caxis([-1.8 -0.4]);
c.Label.String = 'Temperature [{\circ}C]';c.Label.FontSize = 18;
cmocean('thermal');
xlabel('Longitude');ylabel('Depth [m]');
hold on;
[ccont,hcont]=contour(ll_custom_p,-z_p(ii_p(is_p),:),dens_c_annual-1000,cont_lvl_up,'ShowText','on','LineColor','k','LineWidth',0.2);
hcont.LevelList=round(hcont.LevelList,1);
clabel(ccont,hcont,'labelspacing',10.0e+100,'FontWeight','bold','Color','w','FontSize',15);
[ccont,hcont]=contour(ll_custom_p,-z_p(ii_p(is_p),:),dens_c_annual-1000,cont_lvl_down,'ShowText','on','LineColor','k','LineWidth',0.2);
hcont.LevelList=round(hcont.LevelList,1);
clabel(ccont,hcont,'labelspacing',10.0e+10,'FontWeight','bold','Color','w','FontSize',15);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',20,'FontWeight','Bold', 'LineWidth', 2);


ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(2);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),s_annual_pfa);%choose exp here- exp a/b 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;
caxis([33 34]);
c.Label.String = 'Salinity [psu]';c.Label.FontSize = 28;
cmocean('haline');
xlabel('Longitude');ylabel('Depth [m]');
hold on;
[ccont,hcont]=contour(ll_custom_p,-z_p(ii_p(is_p),:),dens_c_annual-1000,cont_lvl_up,'ShowText','on','LineColor','k','LineWidth',0.2);
hcont.LevelList=round(hcont.LevelList,1);
clabel(ccont,hcont,'labelspacing',10.0e+100,'FontWeight','bold','Color','r','FontSize',15);
[ccont,hcont]=contour(ll_custom_p,-z_p(ii_p(is_p),:),dens_c_annual-1000,cont_lvl_down,'ShowText','on','LineColor','k','LineWidth',0.2);
hcont.LevelList=round(hcont.LevelList,1);
clabel(ccont,hcont,'labelspacing',10.0e+100,'FontWeight','bold','Color','r','FontSize',15);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',20,'FontWeight','Bold', 'LineWidth', 2);


%% 5b. Winter - relative to annual - Experiment A/B

%% Temperature

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(3);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),t_winter_pfa-t_annual_pfa); %choose exp here- exp a/b
%hold on;plot(lon(ii_p(is_p)),-h(ii_p(is_p)));plot(lon(ii_p(is_p)),-zi(ii_p(is_p))); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;
caxis([-0.2 0.2]);
c.Label.String = 'Temperature [{\circ}C]';c.Label.FontSize = 28;
%cmocean('thermal');
cmocean('balance','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),t_winter_pfa-t_annual_pfa,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);

%% Salinity

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(4);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),s_winter_pfa-s_annual_pfa);%choose exp here- exp a/b
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;
caxis([-0.15 0.15]);
c.Label.String = 'Salinity [psu]';c.Label.FontSize = 28;
%cmocean('haline');
cmocean('delta','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),s_winter_pfa-s_annual_pfa,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);


%% 5c. Summer - relative to annual - Experiment A/B

%% Temperature 

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(5);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),t_summer_pfa-t_annual_pfa);%choose exp here- exp a/b
%hold on;plot(lon(ii_p(is_p)),-h(ii_p(is_p)));plot(lon(ii_p(is_p)),-zi(ii_p(is_p))); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;
caxis([-0.2 0.2]);
c.Label.String = 'Temperature [{\circ}C]';c.Label.FontSize = 28;
%cmocean('thermal');
cmocean('balance','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),t_summer_pfa-t_annual_pfa,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);

%% Salinity

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(6);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),s_summer_pfa-s_annual_pfa);%choose exp here- exp a/b
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;
caxis([-0.15 0.15]);
c.Label.String = 'Salinity [psu]';c.Label.FontSize = 28;
%cmocean('haline');
cmocean('delta','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),s_summer_pfa-s_annual_pfa,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);


%% 5d. Winter - Relative to standard run - exp A

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(7);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),t_winter_pfa-t_winter_pf);
%hold on;plot(lon(ii_p(is_p)),-h(ii_p(is_p)),'k','LineWidth',1);plot(lon(ii_p(is_p)),-zi(ii_p(is_p)),'k','LineWidth',1); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;caxis([-0.1 0.1]);
c.Label.String = 'Temperature [{\circ}C]';c.Label.FontSize = 28;
cmocean('balance','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),t_winter_pfa-t_winter_pf,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(8);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),s_winter_pfa-s_winter_pf); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;caxis([-0.1 0.1]);colormap jet;
c.Label.String = 'Salinity [psu]';c.Label.FontSize = 28;
cmocean('delta','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),s_winter_pfa-s_winter_pf,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);


%% 5e. Summer - Relative to standard run - exp A

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(9);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),t_summer_pfa-t_summer_pf);
%hold on;plot(lon(ii_p(is_p)),-h(ii_p(is_p)),'k','LineWidth',1);plot(lon(ii_p(is_p)),-zi(ii_p(is_p)),'k','LineWidth',1); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;caxis([-0.1 0.1]);colormap jet;
c.Label.String = 'Temperature [{\circ}C]';c.Label.FontSize = 28;
cmocean('balance','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),t_summer_pfa-t_summer_pf,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);


ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(10);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),s_summer_pfa-s_summer_pf);
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;caxis([-0.1 0.1]);colormap jet;
c.Label.String = 'Salinity [psu]';c.Label.FontSize = 28;
cmocean('delta','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),s_summer_pfa-s_summer_pf,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);


%% 5f. Winter - exp B - exp A

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(11);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),t_winter_pfb-t_winter_pfa);
%hold on;plot(lon(ii_p(is_p)),-h(ii_p(is_p)),'k','LineWidth',1);plot(lon(ii_p(is_p)),-zi(ii_p(is_p)),'k','LineWidth',1); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;caxis([-0.1 0.1]);
c.Label.String = 'Temperature [{\circ}C]';c.Label.FontSize = 28;
cmocean('balance','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),t_winter_pfb-t_winter_pfa,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(12);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),s_winter_pfb-s_winter_pfa); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;caxis([-0.1 0.1]);colormap jet;
c.Label.String = 'Salinity [psu]';c.Label.FontSize = 28;
cmocean('delta','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),s_winter_pfb-s_winter_pfa,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);


%% 5g. Summer - exp B - exp A

ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(13);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),t_summer_pfb-t_summer_pfa);
%hold on;plot(lon(ii_p(is_p)),-h(ii_p(is_p)),'k','LineWidth',1);plot(lon(ii_p(is_p)),-zi(ii_p(is_p)),'k','LineWidth',1); 
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;caxis([-0.1 0.1]);colormap jet;
c.Label.String = 'Temperature [{\circ}C]';c.Label.FontSize = 28;
cmocean('balance','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),t_summer_pfb-t_summer_pfa,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);


ll_custom_p = repmat(lon(ii_p(is_p)),[1 23]);
figure(14);clf
pcolorjw(ll_custom_p,-z_p(ii_p(is_p),:),s_summer_pfb-s_summer_pfa);
set(gca,'Color',[0.7 0.7 0.7]);c=colorbar;caxis([-0.1 0.1]);colormap jet;
c.Label.String = 'Salinity [psu]';c.Label.FontSize = 28;
cmocean('delta','pivot',0);
xlabel('Longitude');ylabel('Depth [m]');
hold on;
contour(ll_custom_p,-z_p(ii_p(is_p),:),s_summer_pfb-s_summer_pfa,[0 0],'ShowText','off','LineColor','k','LineWidth',0.5);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',28,'FontWeight','Bold', 'LineWidth', 2);


