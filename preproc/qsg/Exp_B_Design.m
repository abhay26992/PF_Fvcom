%% Petermann Fjord experiment B design

% Using Exp_B_Design_Plan_btrials, we have selected the <b8> configuration
% to represent experiment B.

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;

%% Section 1 - East - West cross section along (near) the GL

load('IBN_M.mat')
figure(1);
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;caxis([0 500]);hold on

x1=-60.45+360;
y1=80.52;
x2=-59.59+360;
y2=80.62;

xn = 2000;
xx = linspace(x1,x2,xn);
yy = linspace(y1,y2,xn);

for i = 1:numel(xx)
    dum = spheredist(yy(i),xx(i)-360,Mobj.lat,Mobj.lon-360);
    is(i) = find(dum == min(dum));
end
is = unique(is,'stable');

scatter(Mobj.lon(is),Mobj.lat(is),100,'k','o','filled');hold off;

h_cross1=Mobj.h(is);h_crosslat_1=Mobj.lat(is);h_crosslon_1=Mobj.lon(is);

% For cross-section
zcrs = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-Mobj.h),Mobj.siglay),Mobj.zisf); 
zcrs = zcrs(is,:);
ll = repmat(h_crosslon_1,[1 23]);
dval=ones(size(ll));

figure(100);clf
pcolorjw(ll,-zcrs,dval);
set(gca,'Color','k');
ylim([-1000 0]);

%% Section 2 - East - West cross section (further downstream) along the GL

load('IBN_M.mat')
figure(2);
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;caxis([0 500]);hold on

x3=-60.78+360;
y3=80.58;
x4=-59.95+360;
y4=80.69;

xn2 = 2000;
xx2 = linspace(x3,x4,xn2);
yy2 = linspace(y3,y4,xn2);

for i = 1:numel(xx2)
    dum = spheredist(yy2(i),xx2(i)-360,Mobj.lat,Mobj.lon-360);
    is2(i) = find(dum == min(dum));
end
is2 = unique(is2,'stable');

scatter(Mobj.lon(is2),Mobj.lat(is2),100,'k','o','filled');hold off;

h_cross2=Mobj.h(is2);h_crosslat_2=Mobj.lat(is2);h_crosslon_2=Mobj.lon(is2);

% For cross-section
zcrs2 = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-Mobj.h),Mobj.siglay),Mobj.zisf); %to change to old_h : use Mobj.h instead of new_h
zcrs2 = zcrs2(is2,:);
ll2 = repmat(h_crosslon_2,[1 23]);
dval2=ones(size(ll2));

figure(200);clf
pcolorjw(ll2,-zcrs2,dval2);
set(gca,'Color','k');
ylim([-1000 0]);


%% Line Plot

figure(3);clf
plot((Mobj.zisf(is)).*(-1),'c');
hold on;plot((h_cross1).*(-1),'r');
plot((Mobj.zisf(is2)).*(-1),'b');
plot((h_cross2).*(-1),'m');
grid on;ylim([-1000 0]);xlim([0 126]);
xlabel('x');ylabel('z [m]');

%% Line Plot for @GL nodes

load ./xr; load ./yr;
ind_gl=find(inpolygon(Mobj.x,Mobj.y,xr,yr));
h_gl=Mobj.h(ind_gl);z_gl=Mobj.zisf(ind_gl);

figure(4);clf
plot(1:numel(ind_gl),flipud(z_gl.*(-1)),'b');hold on;
plot(1:numel(ind_gl),flipud(h_gl.*(-1)),'r');
grid on;ylim([-1000 0]);xlim([0 89]);
xlabel('x');ylabel('z [m]');

%% Extract cell numbers of interest

%See Exp_B_Design_Plan_btrials.m for details.

%ch_id=[25;36;44;55;67;80;92;98;112]; %exp. id - b6 
ch_id=[18;36;44;55;67;80;92;98;112]; %exp. id - b8 
%ch_id=[36;44;55;67;80;92;98;112]; %exp. id - b9 

is_ch=is(ch_id);%is nodes @channel id [can do the same using is2 instead]

%Visualise on zisf map
figure(5);clf
plot_field(Mobj,Mobj.zisf);
colormap jet;caxis([0 500]);hold on;
scatter(Mobj.x(is_ch),Mobj.y(is_ch),100,'k','o','filled');
%Add the GL river cells (glc_bcell) on the plot to select the corresponding
%cells at the GL (simply eyeball it for the time being). 
load('./tge.mat')
bcell=find(isbce==1);%find all cells that are on a solid boundary (not obc)
glc_1=find(ismember(Mobj.tri(:,1),ind_gl));
glc_2=find(ismember(Mobj.tri(:,2),ind_gl));
glc_3=find(ismember(Mobj.tri(:,3),ind_gl));
glc=sort([glc_1;glc_2;glc_3]);
glc_bcell_ind=find(ismember(glc,bcell));
glc_bcell=unique(glc(glc_bcell_ind));
scatter(Mobj.xc(glc_bcell),Mobj.yc(glc_bcell),100,'c','o','filled');

%% Make the selection

%glr_x=[-2.7e+05;-2.685e+05;-2.676e+05;-2.663e+05;-2.648e+05;-2.628e+05;-2.61e+05;-2.598e+05;-2.574e+05]; %b6
%glr_y=[-9.923e+05;-9.914e+05;-9.909e+05;-9.903e+05;-9.897e+05;-9.894e+05;-9.891e+05;-9.888e+05;-9.878e+05]; %b6

glr_x=[-2.707e+05;-2.685e+05;-2.676e+05;-2.663e+05;-2.648e+05;-2.628e+05;-2.61e+05;-2.598e+05;-2.574e+05]; %b8
glr_y=[-9.927e+05;-9.914e+05;-9.909e+05;-9.903e+05;-9.897e+05;-9.894e+05;-9.891e+05;-9.888e+05;-9.878e+05]; %b8

%glr_x=[-2.685e+05;-2.676e+05;-2.663e+05;-2.648e+05;-2.628e+05;-2.61e+05;-2.598e+05;-2.574e+05]; %b9
%glr_y=[-9.914e+05;-9.909e+05;-9.903e+05;-9.897e+05;-9.894e+05;-9.891e+05;-9.888e+05;-9.878e+05]; %b9

%% Visualise
scatter(glr_x,glr_y,100,'b','o','filled');
%Get the cells from spheredist
[glr_lat,glr_lon]=polarstereo_inv(glr_x,glr_y);
[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);

for i = 1:numel(glr_lon)
    dum = spheredist(glr_lat(i),glr_lon(i)-360,latc,lonc-360);
    rc_expb(i) = find(dum == min(dum));
end
rc_expb = (unique(rc_expb,'stable'))';

%Visualise again..and save..
scatter(Mobj.xc(rc_expb),Mobj.yc(rc_expb),100,'m','o','filled');

save rc_expb rc_expb; %save the 9 selected river cells for exp B


