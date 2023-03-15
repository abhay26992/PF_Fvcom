%% Up-sample the nesting variables

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));

%% Get/Create the hourly time stamp 

%%

%fn='pf_n_to_jul_el_obc.nc'; %Using the merged (A4+AOTIM) tidal file
%hourly=ncread(fn,'time'); 

% Added to create a vector of hourly time-series for the large FVCOM runs:

% year=2015;
% month_in=01;month_fn=12;
% day_in=01;day_fn=31;
% 
% start=datenum(year,month_in,day_in);
% final=datenum(year,month_fn,day_fn);
% ts=datenum(year,month_in,day_in,0:(final-start)*24.',0,0)';
% check_ts=datestr(ts);

%% July 01, 2014 , 00:00:00 - Jan 01, 2015, 00:00:00

% year_in=2014;year_fn=2015;
% month_in=07;month_fn=01;
% day_in_fn=01;
% 
% start=datenum(year_in,month_in,day_in_fn);
% final=datenum(year_fn,month_fn,day_in_fn);
% ts=datenum(year_in,month_in,day_in_fn,0:(final-start)*24.',0,0)';
% check_ts=datestr(ts);

%% Jan 01, YYYY , 00:00:00 - Jan 01, YYYY+1, 00:00:00

year_in=2015;year_fn=2016;
month_in_fn=01;
day_in_fn=01;

start=datenum(year_in,month_in_fn,day_in_fn);
final=datenum(year_fn,month_in_fn,day_in_fn);
ts=datenum(year_in,month_in_fn,day_in_fn,0:(final-start)*24.',0,0)';
check_ts=datestr(ts);

%%
%Save the time-series to make tidal predictions using AOTIM on Stallo:
%save ts_jul_dec_2014 ts 

%Convert to FVCOM time:
hourly=ts-678942;
check_hourly=datestr(double(hourly));
%save jan012017_jan012018.mat hourly

%% Get the daily time stamp and all the nested variables

fn1='pf_jan_dec_15.nc';
%ncdisp(fn1)
daily=ncread(fn1,'time');

% Nesting variables
zeta=ncread(fn1,'zeta');
temp=ncread(fn1,'temp');
salinity=ncread(fn1,'salinity');
u=ncread(fn1,'u');v=ncread(fn1,'v');
ua=ncread(fn1,'ua');va=ncread(fn1,'va');

%% Size : Node/Nele x Time
hs_zeta=zeros(length(zeta),length(hourly));
hs_ua=zeros(length(ua),length(hourly));
hs_va=zeros(length(va),length(hourly));

for p=1:size(zeta,1)
    hs_zeta(p,:)=interp1(daily',zeta(p,:),hourly');
end

% Add tidal elevations from the AOTIM 2018 tidal solution:
load nesting_tidal_elevation_jan_dec_15_4const.mat
hs_zeta=hs_zeta+zt;

for p=1:size(ua,1)
    hs_ua(p,:)=interp1(daily',ua(p,:),hourly');
    hs_va(p,:)=interp1(daily',va(p,:),hourly');
end

%% Size : Node/Nele x Siglay x Time
hs_temp=zeros(length(temp),size(temp,2),length(hourly));
hs_salinity=zeros(length(salinity),size(salinity,2),length(hourly));

for p=1:size(temp,1)
    for q=1:size(temp,2)
        hs_temp(p,q,:)=interp1(daily',squeeze(temp(p,q,:)),hourly');
        hs_salinity(p,q,:)=interp1(daily',squeeze(salinity(p,q,:)),hourly');
    end
end

hs_u=zeros(length(u),size(u,2),length(hourly));
hs_v=zeros(length(v),size(v,2),length(hourly));

for p=1:size(u,1)
    for q=1:size(u,2)
        hs_u(p,q,:)=interp1(daily',squeeze(u(p,q,:)),hourly');
        hs_v(p,q,:)=interp1(daily',squeeze(v(p,q,:)),hourly');
    end
end

%% Other variables:

hyw=ncread(fn1,'hyw');
weight_cell=ncread(fn1,'weight_cell');
weight_node=ncread(fn1,'weight_node');


hs_hyw=zeros(length(hyw),size(hyw,2),length(hourly));

for p=1:size(hyw,1)
    for q=1:size(hyw,2)
        hs_hyw(p,q,:)=interp1(daily',squeeze(hyw(p,q,:)),hourly');
    end
end


hs_weight_cell=zeros(length(weight_cell),length(hourly));
hs_weight_node=zeros(length(weight_node),length(hourly));


for p=1:size(weight_cell,1)
    hs_weight_cell(p,:)=interp1(daily',weight_cell(p,:),hourly');
end

for p=1:size(weight_node,1)
    hs_weight_node(p,:)=interp1(daily',weight_node(p,:),hourly');
end


% Nesting time/Itime/Itime2:

%hs_time=ncread(fn,'time');hs_Itime=ncread(fn,'Itime');hs_Itime2=ncread(fn,'Itime2');

% Read hs_time from the externally created vector of hourly time steps:
hs_time=hourly;
%% All unchanged variables (same as in the pf_new_topo_nest.nc)

nv=ncread(fn1,'nv');
x=ncread(fn1,'x');y=ncread(fn1,'y');
xc=ncread(fn1,'xc');yc=ncread(fn1,'yc');
lon=ncread(fn1,'lon');lat=ncread(fn1,'lat');
lonc=ncread(fn1,'lonc');latc=ncread(fn1,'latc');
h=ncread(fn1,'h');h_center=ncread(fn1,'h_center');
siglay=ncread(fn1,'siglay');siglev=ncread(fn1,'siglev');
siglay_center=ncread(fn1,'siglay_center');siglev_center=ncread(fn1,'siglev_center');

%% Save all variables to a new struct -> (hf_nest_st)

hf_nest_st.nv=nv;
hf_nest_st.x=x;hf_nest_st.y=y;
hf_nest_st.xc=xc;hf_nest_st.yc=yc;
hf_nest_st.lon=lon;hf_nest_st.lat=lat;
hf_nest_st.lonc=lonc;hf_nest_st.latc=latc;
hf_nest_st.h=h;
hf_nest_st.siglay=siglay;hf_nest_st.siglev=siglev;


hf_nest_st.hs_time=hs_time;
hf_nest_st.hs_salinity=hs_salinity;
hf_nest_st.hs_temp=hs_temp;
hf_nest_st.hs_hyw=hs_hyw;
hf_nest_st.hs_zeta=hs_zeta;
hf_nest_st.hs_u=hs_u;
hf_nest_st.hs_v=hs_v;
hf_nest_st.hs_ua=hs_ua;
hf_nest_st.hs_va=hs_va;
hf_nest_st.hs_weight_cell=hs_weight_cell;
hf_nest_st.hs_weight_node=hs_weight_node;


%% Write all variables to a netcdf file (The New Nesting File):

load ./ngrd.mat; load ./M.mat; nest_struct=hf_nest_st;
fileprefix='pf_hf_jan_dec_15';

write_hf_nest(ngrd,Mobj,nest_struct,fileprefix);

% View the new netcdf file:
ncdisp([fileprefix, '.nc'])

%% Plot:

% From the hourly AOTIM+A4 jul_el_obc.nc:
%elevation=ncread(fn,'elevation');
%hourly=ncread(fn,'time'); %in FVCOM time
%hourly=hourly+678942; %Add time-offset to go back to Matlab convention


% From the new nesting file (hourly):
elev_nest=ncread([fileprefix, '.nc'],'zeta');
nest_time=ncread([fileprefix, '.nc'],'time'); %in FVCOM time
nest_time=nest_time+678942; %Add time-offset to go back to Matlab convention

% From the old nesting file (daily):
daily=ncread(fn1,'time');
daily=daily+678942; %Add time-offset to go back to Matlab convention
zeta=ncread(fn1,'zeta');


figure(1);clf
%plot(hourly,elevation(100,:),'linewidth',2);hold on;
%scatter(nest_time,elev_nest(100,:),25,'r');hold on;
%scatter(daily,zeta(100,:),150,'k','filled');
plot(nest_time,elev_nest(100,:),'linewidth',2);hold on;
scatter(daily,zeta(100,:),150,'k','filled');
datetick('x','ddmmyyyy','keeplimits');
xlabel('Time');ylabel('Elevation');
%legend('jul el obc.nc','old nest.nc','new nest.nc');
legend('new nest.nc','old nest.nc');
title('Elevation time series - Node 100')


% figure(2);clf
% plot(hourly,elevation(100,:),'b','linewidth',2);hold on;
% plot(nest_time,elev_nest(100,:),'r','linewidth',2);
% datetick('x','ddmmyyyy','keeplimits');
% xlabel('Time');ylabel('Elevation');
% legend('jul el obc.nc','new nest.nc');
% title('Elevation time series - Node 100')

figure(3);clf
t=ncread([fileprefix, '.nc'],'t');checknan=find(isnan(t)==1);
t=double(t);t1=squeeze(t(100,1,:));
plot(ts(end-48:end),t1(end-48:end),'linewidth',2);grid on;datetick('x',31);
