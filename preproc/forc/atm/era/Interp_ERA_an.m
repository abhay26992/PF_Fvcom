clear all;close all;clc;

%%

year=2016;

load ./erawdind  % only uses the era data within the model domain
use erawdind
irl=length(idr);
jrl=length(jdr);


% interpolating index and coefficients
load ./N4ERA
use N4ERA

nn=length(nindex);
nc=length(cindex);

%% For Figures:

lone=ncread('era_pmf_an_2010_19.nc','longitude');late=ncread('era_pmf_an_2010_19.nc','latitude');
late=late';
[lat,lon]=meshgrid(late,lone);
lon=lon(idr,jdr);lat=lat(idr,jdr);

%% Read in ERA data, pay attention to the signs and units%%%%%%%%%%%

%% 2 - Meter Air Temperature -> Surface air temperature in Degree C

folder='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/';
file='era_pmf_an_2010_19.nc';


% Change time for ERA (-Hours- since "1900-01-01")
era_time=ncread([folder file],'time'); %hours since 1900-01-01 00:00:00
T = datetime(1900,1,1) + hours(era_time);
time_era=datenum(T);

% Set new time-limit
time_lim=[datenum(2016,01,01) datenum(2017,01,01)];
ind1=find(time_era>=time_lim(1));
ind2=find(time_era<=time_lim(2));
ind=intersect(ind1,ind2);
time_era=time_era(ind);check_time=datestr(time_era);
nt=length(time_era);


% Dimension is different for ERA (compared to RACMO)
T2m=ncread([folder file],'t2m',[idr(1) jdr(1) 1],[irl jrl Inf]);
T2m=T2m-273.15; % Convert to Celsius
T2m=T2m(:,:,ind); %Dimension is different for ERA (compared to RACMO)


for i=1:nt
      
    TT=squeeze(T2m(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    TA_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

%Plot to check:
%grid_north_pole_longitude = 142.5;grid_north_pole_latitude  = 18;
%SP_Lon=grid_north_pole_longitude-180;SP_Lat=grid_north_pole_latitude*(-1);
%SP_coor=[SP_Lon SP_Lat];
%rlon=ncread('t2m.nc','rlon');rlat=ncread('t2m.nc','rlat');

figure(1);clf
%Plot cropped long-lat and data
pcolor(lon,lat,T2m(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,100,TA_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/TA_ERA_' num2str(year)],'time_era','TA_ERA','-v7.3')
clear TA_ERA T2m

%% Surface Air Pressure

Psurf=ncread([folder file],'sp',[idr(1) jdr(1) 1],[irl jrl Inf]);
Psurf=Psurf(:,:,ind); %Dimension is different for ERA (compared to RACMO)


for i=1:nt
      
    TT=squeeze(Psurf(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    PA_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

figure(2);clf
%Plot cropped long-lat and data
pcolor(lon,lat,Psurf(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,100,PA_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/PA_ERA_' num2str(year)],'time_era','PA_ERA','-v7.3')
clear PA_ERA Psurf

%% U-wind component

Uwind=ncread([folder file],'u10',[idr(1) jdr(1) 1],[irl jrl Inf]);
Uwind=Uwind(:,:,ind); %Dimension is different for ERA (compared to RACMO)

UW_ERA=NaN(nc,nt);

for i=1:nt
      
    TT=squeeze(Uwind(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    UW_ERA(:,i)=sum(TT(cindex).*ccoef,2);
    clear TT
end 

figure(3);clf
%Plot cropped long-lat and data
pcolor(lon,lat,Uwind(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M;[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
scatter(lonc,latc,100,UW_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/UW_ERA_' num2str(year)],'time_era','UW_ERA','-v7.3')
clear UW_ERA Uwind

%% V-wind component

Vwind=ncread([folder file],'v10',[idr(1) jdr(1) 1],[irl jrl Inf]);
Vwind=Vwind(:,:,ind); %Dimension is different for ERA (compared to RACMO)

VW_ERA=NaN(nc,nt);

for i=1:nt
      
    TT=squeeze(Vwind(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    VW_ERA(:,i)=sum(TT(cindex).*ccoef,2);
    clear TT
end 

figure(4);clf
%Plot cropped long-lat and data
pcolor(lon,lat,Vwind(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M;[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
scatter(lonc,latc,100,VW_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/VW_ERA_' num2str(year)],'time_era','VW_ERA','-v7.3')
clear VW_ERA Vwind

%% Specific and Relative Humidity:

%% Specific Humidity:

% Use dew-point temperature and surface pressure to calculate specific
% humidity:

% Get the dew point temperature:
t_dew=ncread([folder file],'d2m',[idr(1) jdr(1) 1],[irl jrl Inf]);
t_dew=t_dew(:,:,ind); %Dimension is different for ERA (compared to RACMO)
t_dew=t_dew-273.15; %Convert to Degree Celsius

% Get the surface pressure:
Psurf=ncread([folder file],'sp',[idr(1) jdr(1) 1],[irl jrl Inf]);
Psurf=Psurf(:,:,ind); %Dimension is different for ERA (compared to RACMO)
Psurf=Psurf./100; %Convert Pascal to mbar

% Calculate the saturation water vapour pressure (Tetens formula)

e = 6.112*exp((17.67*t_dew)./(t_dew + 243.5)); %in mbar

% Calculate the specific humidity (using dew point temperature and
% saturation water vapour pressure):
q = (0.622 * e)./(Psurf - (0.378 * e)); % q in g/kg

%%%%%%% Convert to kg/kg (for A4-CICE forcing)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_k=q./1000;

for i=1:nt
      
    TT=squeeze(q_k(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    SH_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

figure(5);clf
%Plot cropped long-lat and data
pcolor(lon,lat,q_k(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,100,SH_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/SH_ERA_' num2str(year)],'time_era','SH_ERA','-v7.3')
clear SH_ERA q_k

%% Relative Humidity:

% Get the 2-m air temperature 
T2m=ncread([folder file],'t2m',[idr(1) jdr(1) 1],[irl jrl Inf]);
T2m=T2m(:,:,ind); %Dimension is different for ERA (compared to RACMO)
T2m=T2m-273.15; % Convert to Celsius


es = 6.112*exp((17.67*T2m)./(T2m + 243.5)); % Saturation Vapor Pressure in mb
e = 6.112*exp((17.67*t_dew)./(t_dew + 243.5)); % Vapour Pressure in mb

%RH = 100.0 * (e./es); % Relative Humidity in percent
RH= e./es;

for i=1:nt
      
    TT=squeeze(RH(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    RH_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

figure(6);clf
%Plot cropped long-lat and data
pcolor(lon,lat,RH(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,100,RH_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/RH_ERA_' num2str(year)],'time_era','RH_ERA','-v7.3')
clear RH_ERA RH

%% Cloud - Cover (for A4-CICE)

tcc=ncread([folder file],'tcc',[idr(1) jdr(1) 1],[irl jrl Inf]);
tcc=tcc(:,:,ind); %Dimension is different for ERA (compared to RACMO)

for i=1:nt
      
    TT=squeeze(tcc(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    TCC_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

figure(7);clf
%Plot cropped long-lat and data
pcolor(lon,lat,tcc(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,100,TCC_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/TCC_ERA_' num2str(year)],'time_era','TCC_ERA','-v7.3')
clear TCC_ERA tcc