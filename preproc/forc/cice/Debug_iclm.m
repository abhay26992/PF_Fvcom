
clear all;close all;clc;

%% FVCOM

load ./M.mat

aice=ncread('pf_jan_dec_15_icenudge.nc','AICE');aice=double(aice);
hice=ncread('pf_jan_dec_15_icenudge.nc','HICE');hice=double(hice);

uuice=ncread('pf_jan_dec_15_icenudge.nc','UICE');uuice=double(uuice);
vvice=ncread('pf_jan_dec_15_icenudge.nc','VICE');vvice=double(vvice);

%% A4 - Grid
load rind
use rind

irl=length(idr);
jrl=length(jdr);

iul=length(idu);
jul=length(jdu);

ivl=length(idv);
jvl=length(jdv);

a4lon=ncread('iceh.2014-07-01.nc','TLON',[idr(1) jdr(1)],[irl jrl]);
a4lat=ncread('iceh.2014-07-01.nc','TLAT',[idr(1) jdr(1)],[irl jrl]);
%a4mask=ncread('iceh.2014-07-01.nc','blkmask',[idr(1) jdr(1)],[irl jrl]);
a4vice=ncread('iceh.2014-07-01.nc','vicen_d',[idr(1) jdr(1) 1 1],[irl jrl Inf Inf]);a4vice=a4vice(:,:,1);
a4aice=ncread('iceh.2014-07-01.nc','aicen_d',[idr(1) jdr(1) 1 1],[irl jrl Inf Inf]);a4aice=a4aice(:,:,1);

% [c,h]=contour(a4vice);
% idx = c(1,:)<1 ; 
% c(:,idx) = NaN;

%% Plot:

figure(1);clf
plot_field(Mobj,hice(:,1),'coordinate','spherical');colorbar;colormap jet;hold on;
pcolor(a4lon-360,a4lat,a4vice);shading interp;hold on;
%plot(c(1,:),c(2,:),'r')

% k=boundary(Mobj.lon,Mobj.lat,1); 
% b_fvlon=Mobj.lon(k);b_fvlat=Mobj.lat(k);
% figure(2);clf
% plot(b_fvlon,b_fvlat,'k','LineWidth',2);hold on;

figure(2);clf
plot_field(Mobj,aice(:,300),'coordinate','spherical');colorbar;colormap jet;hold on; %tstep = 60,180,300

figure(3);clf
plot_field(Mobj,hice(:,300),'coordinate','spherical');colorbar;colormap jet;hold on; %tstep = 60,180,300

cv=uuice+1i*vvice;
figure(4);clf
spd=elems2nodes(abs(cv(:,300)),Mobj); %tstep = 60,180,300
plot_field(Mobj,spd,'coordinate','spherical');colormap jet;caxis([0 1.0]);