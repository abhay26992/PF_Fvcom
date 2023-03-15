%% Write out the A4 (climatology) nesting file

% Tides are to be kept the same (Upsampled daily A4 + hourly AOTIM). For
% other variables, we use monthly climatologies. As a result, there will be 
% two different time-steps, i.e. 12+1 monthly time-steps and ~365*24 hourly 
% time-steps in the nesting file. The easiest way to implement this is to 
% interpolate the monthly climatology-steps on to the hourly-steps.

% Split in 3 parts. This is part (3) of (3). 

% See previous --> A4_clima.m (Part 1) and A4_merge_clima.m (Part 2)

% Written by Abhay Prakash (abhay.prakash@natgeo.su.se)

clear all;close all;clc;
dir=pwd;
load ./A4_2007_2017/Boundary/2007-09/clima_0708;

%% Hourly Tides:

% There are several ways to implement this. Loading directly from
% previously generated nesting files is one approach. However: 

% For 2017, A4 variables were available only until December 01. With the
% climatology approach, this can be extended till January 01. We need to
% revisit the AOTIM data to generate the 4-constituent tidal solution from
% Jan 01 2017 - Jan 01 2018, so as to have 3 full years of simulation data. 

% In the approach below, we load the pregenerated hourly time steps and 
% interpolate the A4 daily zeta to the hourly steps. The pregenerated AOTIM
% tides are loaded and added to the upsampled zeta. 

cd C:\PhD\FVCOM\Matlab_Repository\Petermann_Bathy\Setup_Nesting\ssh_adjustment\Nesting_Files\AOTIM_Time_Series\V2.0
load ./jan012014_jan012015.mat;
hourly=hourly+678942; %Convert from FVCOM time


cd C:\PhD\FVCOM\Matlab_Repository\Petermann_Bathy\Setup_Nesting\ssh_adjustment\Nesting_Files\Nesting_Daily\V2.0
fn='pf_jan_dec_14.nc';
daily=ncread(fn,'time');
daily=daily+678942; %Convert from FVCOM time

zeta=ncread(fn,'zeta');
hs_zeta=zeros(length(zeta),length(hourly));
for p=1:size(zeta,1)
    hs_zeta(p,:)=interp1(daily',zeta(p,:),hourly');
end

% Add tidal elevations from the AOTIM 2018 tidal solution:
cd C:\PhD\FVCOM\Matlab_Repository\Petermann_Bathy\Setup_Nesting\ssh_adjustment\Nesting_Files\AOTIM_Time_Series\V2.0
load ./nesting_tidal_elevation_jan_dec_14_4const.mat
hs_zeta=hs_zeta+zt;


cd(dir);
%% Interpolate monthly climatology variables to hourly time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Append the January climatology to the end of the struct 

% The simulation runs from Jan-01-YYYY to Jan-01-YYYY+1 , and therefore, we
% need to copy and paste the January climatology at the end.

%cat will place B after A
clima_0708.temp=cat(3,clima_0708.temp,clima_0708.temp(:,:,1)); 
clima_0708.salt=cat(3,clima_0708.salt,clima_0708.salt(:,:,1));
clima_0708.u=cat(3,clima_0708.u,clima_0708.u(:,:,1));
clima_0708.v=cat(3,clima_0708.v,clima_0708.v(:,:,1));
clima_0708.ubar=cat(2,clima_0708.ubar,clima_0708.ubar(:,1));
clima_0708.vbar=cat(2,clima_0708.vbar,clima_0708.vbar(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. Get the monthly climatology time-steps

year='2014'; % The nesting file year
months=1:12;

for i=1:length(months)
    t_id(i,1)=datenum(str2double(year),months(i),01,0,0,0);
end

% January - 01 of YYYY+1
t=datenum(str2double(year)+1,months(1),01,0,0,0);
t_id=[t_id;t];
check=datestr(t_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 3. Interpolate monthly climatology fields (t_id) to hourly:

%% 3a. Size : Node/Nele x Time

hs_ua=zeros(length(clima_0708.ubar),length(hourly));
hs_va=zeros(length(clima_0708.vbar),length(hourly));

for p=1:size(clima_0708.ubar,1)
    hs_ua(p,:)=interp1(t_id',clima_0708.ubar(p,:),hourly');
    hs_va(p,:)=interp1(t_id',clima_0708.vbar(p,:),hourly');
end

%% 3b. Size : Node/Nele x Siglay x Time
hs_temp=zeros(length(clima_0708.temp),size(clima_0708.temp,2),length(hourly));
hs_salinity=zeros(length(clima_0708.salt),size(clima_0708.salt,2),length(hourly));

for p=1:size(clima_0708.temp,1)
    for q=1:size(clima_0708.temp,2)
        hs_temp(p,q,:)=interp1(t_id',squeeze(clima_0708.temp(p,q,:)),hourly');
        hs_salinity(p,q,:)=interp1(t_id',squeeze(clima_0708.salt(p,q,:)),hourly');
    end
end

hs_u=zeros(length(clima_0708.u),size(clima_0708.u,2),length(hourly));
hs_v=zeros(length(clima_0708.v),size(clima_0708.v,2),length(hourly));

for p=1:size(clima_0708.u,1)
    for q=1:size(clima_0708.u,2)
        hs_u(p,q,:)=interp1(t_id',squeeze(clima_0708.u(p,q,:)),hourly');
        hs_v(p,q,:)=interp1(t_id',squeeze(clima_0708.v(p,q,:)),hourly');
    end
end

%% Save as struct

clima.u=hs_u;
clima.v=hs_v;
clima.salt=hs_salinity;
clima.temp=hs_temp;
clima.ubar=hs_ua;
clima.vbar=hs_va;
clima.zeta=hs_zeta;
clima.time=hourly';

save ./A4_2007_2017/Boundary/2007-09/clima_2014 clima -v7.3

%% Write out using write_FVCOM_nest

clear all;close all;clc;
load ./ngrd; load ./A4_2007_2017/Boundary/2007-09/clima_2014; 
time=clima.time; nest_type=3; geo2xy=0;
fileprefix='pf_clima14';

write_FVCOM_nest(ngrd,fileprefix,clima,time,nest_type,geo2xy)

