clear all;close all;clc;

%%

year=2016;

load ./erawdind  % only uses the wrf data within the model domain
use erawdind
irl=length(idr);
jrl=length(jdr);


% interpolating index and coefficients
load ./N4ERA
use N4ERA

nn=length(nindex);
nc=length(cindex);

%%  Read in ERA data, pay attention to the signs and units      %%%%%%%%%%%

% The ERA Lat-Lon grid:
lone=ncread('era_pmf_fc_2010_19.nc','longitude');late=ncread('era_pmf_fc_2010_19.nc','latitude');
late=late';
[lat,lon]=meshgrid(late,lone);
lon=lon(idr,jdr);lat=lat(idr,jdr);

%% Total Precipitation

folder='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/';
file='era_pmf_fc_2010_19.nc';

accum_tp=ncread([folder file],'tp',[idr(1) jdr(1) 1],[irl jrl Inf]); %in m
accum_tp=accum_tp./(3*60*60); %in m/s

%% Check (TP = Large Scale Precipitation + Convective Precipitation)

% lsp=ncread([folder file],'lsp',[idr(1) jdr(1) 1],[irl jrl Inf]); %in m
% lsp=lsp./(3*60*60); %in m/s
% 
% cp=ncread([folder file],'cp',[idr(1) jdr(1) 1],[irl jrl Inf]); %in m
% cp=cp./(3*60*60); %in m/s
% 
% sum_tp=lsp+cp;

%%

% Change time for ERA (-Hours- since "1900-01-01")
era_time=ncread([folder file],'time'); %hours since 1900-01-01 00:00:00
T = datetime(1900,1,1) + hours(era_time);
time_era=datenum(T);

% Set new time-limit
%time_lim=[datenum(2016,12,31,21,0,0) datenum(2017,12,31,0,0,0)]; %Forecast -> 2017,01,01,0,0,0 = 2017,01,01,0,0,0 - 2016,12,31,21,0,0,0

% New limit to include all time steps between 31 Dec 00:00:00 - 01 Jan 00:00:00
time_lim=[datenum(2015,12,31,21,0,0) datenum(2017,01,01,0,0,0)];
ind1=find(time_era>=time_lim(1));
ind2=find(time_era<=time_lim(2));
ind=intersect(ind1,ind2);

time_era=time_era(ind);check_time=datestr(time_era);
accum_tp=accum_tp(:,:,ind); %Dimension is different for ERA (compared to RACMO)


%%%%%%%%%%%% Set the "constant" negative rainfall values to 0 %%%%%%%%%%%%%
accum_tp(accum_tp<0)=0; 
check=accum_tp(accum_tp<0); % check=Rain(isnan(Rain)); check=~nnz(Rain);


nt=length(time_era)-1;

tp=NaN(size(lat,1),size(lat,2),nt);
cycle=4; % tp accumulates over 4 time steps (from 03:00 to 12:00 and then from 15:00 to 00:00)
k=(3:cycle:nt)';

for i=2:length(ind)
    in_cycle=find(i == (k(:,1)));
    if ~isempty(in_cycle)
        tp(:,:,k(in_cycle)-1)=accum_tp(:,:,k(in_cycle));
    else
        tp(:,:,i-1)=accum_tp(:,:,i)-accum_tp(:,:,i-1);
    end 
end


for i=1:nt
      
    TT=squeeze(tp(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    TP_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

figure(1);clf
%Plot cropped long-lat and data
pcolor(lon,lat,tp(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,750,TP_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/TP_ERA_' num2str(year)],'time_era','TP_ERA','-v7.3')


%% Evaporation:

accum_evap=ncread([folder file],'e',[idr(1) jdr(1) 1],[irl jrl Inf]); %in m
accum_evap=accum_evap./(3*60*60); %in m/s
accum_evap=accum_evap(:,:,ind); %Dimension is different for ERA (compared to RACMO)

evap=NaN(size(lat,1),size(lat,2),nt);
cycle=4; % evap accumulates over 4 time steps (from 03:00 to 12:00 and then from 15:00 to 00:00)
k=(3:cycle:nt)';

for i=2:length(ind)
    in_cycle=find(i == (k(:,1)));
    if ~isempty(in_cycle)
        evap(:,:,k(in_cycle)-1)=accum_evap(:,:,k(in_cycle));
    else
        evap(:,:,i-1)=accum_evap(:,:,i)-accum_evap(:,:,i-1);
    end 
end

for i=1:nt 
    TT=squeeze(evap(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    Evap_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

figure(2);clf
%Plot cropped long-lat and data
pcolor(lon,lat,evap(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,750,Evap_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/Evap_ERA_' num2str(year)],'time_era','Evap_ERA','-v7.3')

%% Downwelling Longwave Radiation:

accum_dlw=ncread([folder file],'strd',[idr(1) jdr(1) 1],[irl jrl Inf]); %Accumulated radiative flux in Jm-2 

% Change time for ERA (-Hours- since "1900-01-01")
era_time=ncread([folder file],'time'); %hours since 1900-01-01 00:00:00
T = datetime(1900,1,1) + hours(era_time);
time_era=datenum(T);

% Set new time-limit
%time_lim=[datenum(2016,12,31,21,0,0) datenum(2017,12,31,0,0,0)];

% New limit to include all time steps between 31 Dec 00:00:00 - 01 Jan 00:00:00
time_lim=[datenum(2015,12,31,21,0,0) datenum(2017,01,01,0,0,0)];
ind1=find(time_era>=time_lim(1));
ind2=find(time_era<=time_lim(2));
ind=intersect(ind1,ind2);

time_era=time_era(ind);check_time=datestr(time_era);
accum_dlw=accum_dlw(:,:,ind); 

nt=length(time_era)-1;

Lw_down=NaN(size(lat,1),size(lat,2),nt);
cycle=4; % physical fluxes accumulate over 4 time steps (from 03:00 to 12:00 and then from 15:00 to 00:00)
k=(3:cycle:nt)';

for i=2:length(ind)
    in_cycle=find(i == (k(:,1)));
    if ~isempty(in_cycle)
        Lw_down(:,:,k(in_cycle)-1)=accum_dlw(:,:,k(in_cycle));
    else
        Lw_down(:,:,i-1)=accum_dlw(:,:,i)-accum_dlw(:,:,i-1);
    end 
end

Lw_down=Lw_down./(3*3600); %Convert from Jm-2 to Wm-2

% Check:
% maxmax(Lw_down(:,:,1))

for i=1:nt
      
    TT=squeeze(Lw_down(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    LW_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

figure(3);clf
%Plot cropped long-lat and data
pcolor(lon,lat,Lw_down(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,100,LW_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/LW_ERA_' num2str(year)],'time_era','LW_ERA','-v7.3')
clear LW_ERA 

%% Shortwave Radiation:

%%%%%%%%%%%%%%%%%%%%% Net Shortwave Radiation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accum_sw=ncread([folder file],'ssr',[idr(1) jdr(1) 1],[irl jrl Inf]); %Accumulated radiative flux in Jm-2 

% Change time for ERA (-Hours- since "1900-01-01")
era_time=ncread([folder file],'time'); %hours since 1900-01-01 00:00:00
T = datetime(1900,1,1) + hours(era_time);
time_era=datenum(T);

% Set new time-limit
%time_lim=[datenum(2017,01,01,3,0,0) datenum(2017,12,31,3,0,0)]; %Account for the phase-lag seen in the diurnal cycles between ERA and RACMO


% New limit to include all time steps between 31 Dec 00:00:00 - 01 Jan 00:00:00
time_lim=[datenum(2016,01,01,3,0,0) datenum(2017,01,01,3,0,0)];
ind1=find(time_era>=time_lim(1));
ind2=find(time_era<=time_lim(2));
ind=intersect(ind1,ind2);

time_era=time_era(ind);check_time=datestr(time_era);
accum_sw=accum_sw(:,:,ind); 

nt=length(time_era);

Sw_net=NaN(size(lat,1),size(lat,2),nt);
cycle=4; % physical fluxes accumulate over 4 time steps (from 03:00 to 12:00 and then from 15:00 to 00:00)
k=(1:cycle:nt)';

for i=1:length(ind)
    in_cycle=find(i == (k(:,1)));
    if ~isempty(in_cycle)
        Sw_net(:,:,k(in_cycle))=accum_sw(:,:,k(in_cycle));
    else
        Sw_net(:,:,i)=accum_sw(:,:,i)-accum_sw(:,:,i-1);
    end 
end

Sw_net=Sw_net./(3*3600); %Convert from Jm-2 to Wm-2
%Sw_net=Sw_net(:,:,2:end); %Account for the phase-lag seen in the diurnal cycles between ERA and RACMO

for i=1:nt    
    TT=squeeze(Sw_net(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    SW_ERA(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
end 

figure(4);clf
%Plot cropped long-lat and data
pcolor(lon,lat,Sw_net(:,:,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,100,SW_ERA(:,1))
colorbar;colormap jet;
ylim([73 87]);

save(['./Era_Forc/2016/SW_ERA_' num2str(year)],'time_era','SW_ERA','-v7.3')
clear SW_ERA
