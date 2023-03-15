clear all;close all;clc;

year=2016;

load ./wdind  % only uses the racmo data within the model domain
use wdind
irl=length(idr);
jrl=length(jdr);

  
% interpolating index and coefficients  
load ./N4RAC
use N4RAC
  
nn=length(nindex);
nc=length(cindex);
 
  
%% Read in racmo data, pay attention to the signs and units%%%%%%%%%

%% Shortwave Radiation (in W/m^2):

%%%%%%%%%%%%%%%%%%%%%%Shortwave-Down%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/';
file='swsd.nc';

% Change time for RACMO (Days since "1950-01-01")
time=ncread([folder file],'time'); %days since 1950-01-01 00:00:00
time=time+datenum(1950,01,01);
nt=length(time);

d_swrad=ncread([folder file],'swsd',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);

% Set new time-limit
time_lim=[datenum(2016,01,01) datenum(2017,01,01)];
ind1=find(time>=time_lim(1));
ind2=find(time<=time_lim(2));
ind=intersect(ind1,ind2);

time=time(ind);check_time=datestr(time);
d_swrad=d_swrad(:,:,1,ind);
nt=length(time);

for i=1:nt
    TT=squeeze(d_swrad(:,:,i));
    TT=TT(:);
    d_SW(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT 
end

figure(1);clf
lon=ncread('swsd.nc','lon');lat=ncread('swsd.nc','lat');
lon=lon(idr,jdr);lat=lat(idr,jdr);
%Plot cropped long-lat and data
pcolor(lon,lat,d_swrad(:,:,1,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,10,d_SW(:,1))
colorbar;colormap jet;
ylim([73 87]); 

%%%%%%%%%%%%%%%%%%%%%%%% Shortwave - Up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file='swsu.nc';
u_swrad=ncread([folder file],'swsu',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);

u_swrad=u_swrad(:,:,1,ind);

for i=1:nt
    TT=squeeze(u_swrad(:,:,i));
    TT=TT(:);
    u_SW(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end

figure(2);clf
%Plot cropped long-lat and data
pcolor(lon,lat,u_swrad(:,:,1,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,10,u_SW(:,1)) %Note: -ve sign for upwelling rad
colorbar;colormap jet;
ylim([73 87]); 

%%%%%%%%%%%%%%%%%%%%%%%% Net - SW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SW = d_SW+u_SW; clear d_swrad u_swrad

figure(3);clf
load ./M
scatter(Mobj.lon,Mobj.lat,10,SW(:,1))
colorbar;colormap jet;
ylim([73 87]); 

save(['./Racmo_Forc/2016/SW_RAC_' num2str(year)],'time','SW','-v7.3')
clear SW

%% Longwave:
file='lwsd.nc';
d_lwrad=ncread([folder file],'lwsd',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);
d_lwrad=d_lwrad(:,:,1,ind);
 for i=1:nt
   
     TT=squeeze(d_lwrad(:,:,i));
     TT=TT(:);
     d_LW(:,i)=sum(TT(nindex).*ncoef,2);
         
 end
 
figure(4);clf
lon=ncread('lwsd.nc','lon');lat=ncread('lwsd.nc','lat');
lon=lon(idr,jdr);lat=lat(idr,jdr);
%Plot cropped long-lat and data
pcolor(lon,lat,d_lwrad(:,:,1,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,10,d_LW(:,1))
colorbar;colormap jet;
ylim([73 87]); 
 
save(['./Racmo_Forc/2016/LW_RAC_' num2str(year)],'time','d_LW','-v7.3')                                                             
 

clear d_LW d_lwrad
