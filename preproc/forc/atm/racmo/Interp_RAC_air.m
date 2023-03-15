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

%% Read in RACMO data, pay attention to the signs and units%%%%%%%%%%%

%% 2 - Meter Air Temperature -> Surface air temperature in Degree C

folder='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/';
file='t2m.nc';


% Change time for RACMO (Days since "1950-01-01")
time=ncread([folder file],'time'); %days since 1950-01-01 00:00:00
time=time+datenum(1950,01,01);
nt=length(time);


% Dimension is different for RACMO
T2m=ncread([folder file],'t2m',[idr(1) jdr(1) 1 1],[irl jrl  1 Inf]);
T2m=T2m-273.15; % Convert to Celsius

% Set new time-limit
time_lim=[datenum(2016,01,01) datenum(2017,01,01)];
ind1=find(time>=time_lim(1));
ind2=find(time<=time_lim(2));
ind=intersect(ind1,ind2);

time=time(ind);check_time=datestr(time);
T2m=T2m(:,:,1,ind); %Dimension is different for RACMO
nt=length(time);

for i=1:nt
      
    TT=squeeze(T2m(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    TA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

%Plot to check:
%grid_north_pole_longitude = 142.5;grid_north_pole_latitude  = 18;
%SP_Lon=grid_north_pole_longitude-180;SP_Lat=grid_north_pole_latitude*(-1);
%SP_coor=[SP_Lon SP_Lat];
%rlon=ncread('t2m.nc','rlon');rlat=ncread('t2m.nc','rlat');

figure(1);clf
lon=ncread('t2m.nc','lon');lat=ncread('t2m.nc','lat');
lon=lon(idr,jdr);lat=lat(idr,jdr);
%Plot cropped long-lat and data
pcolor(lon,lat,T2m(:,:,1,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,10,TA(:,1))

%% Check the final (merged) forcing product:
%load ./T_2m_forc.mat
%scatter(Mobj.lon,Mobj.lat,10,T_2m_forc(:,1))

%%
colorbar;colormap jet;
ylim([73 87]);
save(['./Racmo_Forc/2016/TA_RAC_' num2str(year)],'time','TA','-v7.3')
clear TA T2m

%% Surface Air Pressure (Already in Pa for RACMO)

file='psurf.nc';
Psurf=ncread([folder file],'psurf',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);
Psurf=Psurf(:,:,1,ind);


for i=1:nt
    TT=squeeze(Psurf(:,:,i));
    TT=TT(:);
    PA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
    
end

figure(2);clf
%Plot cropped long-lat and data
pcolor(lon,lat,Psurf(:,:,1,1));shading flat;hold on;
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,10,PA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Racmo_Forc/2016/PA_RAC_' num2str(year)],'time','PA','-v7.3')

clear PA Psurf

%% Relative Humidity (Dimensionless)

file='rh2m.nc';
RH2m=ncread([folder file],'rh2m',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);
RH2m=RH2m(:,:,1,ind);

for i=1:nt
    
    TT=squeeze(RH2m(:,:,i));
    TT=TT(:);
    RH(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
    
end

figure(3);clf
%Plot cropped long-lat and data
pcolor(lon,lat,RH2m(:,:,1,1));shading flat;hold on;
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,10,RH(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Racmo_Forc/2016/RH_RAC_' num2str(year)],'time','RH','-v7.3')
clear RH RH2m

%% Specific Humidity (Already in kg/kg for RACMO)

% file='q2m.nc';
% Q2m=ncread([folder file],'q2m',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);
% Q2m=Q2m(:,:,1,ind);
% 
% for i=1:nt
%     
%     TT=squeeze(Q2m(:,:,i));
%     TT=TT(:);
%     QA(:,i)=sum(TT(nindex).*ncoef,2);
%     clear TT
%     
% end
% 
% 
% save(['QA' num2str(year)],'time','QA','-v7.3')
% clear QA Q2m
