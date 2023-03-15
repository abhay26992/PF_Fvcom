%% Dumping the Atmospheric forcing files for FVCOM:

clear all;close all;clc;

%% Load common inputs to the write_FVCOM_forcing_var function 
load ./M 

load ./Atm_Forc_T_Vec_2016.mat  % 3-hr time steps from RACMO
time = RACMO_time; clear RACMO_time; 
time=time-datenum(1858,11,17,0,0,0);

fileprefix='mergedforc16';

%% U-V wind velocity:

%%%%%%%%%%%%%%%%%%%%%%%%%%%% U-10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grdtype='cell';
var='uwnd';

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/U10m_forc_16.mat %The Merged Product
data=U10m_forc; clear U10m_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V-10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var='vwnd';

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/V10m_forc_16.mat %The Merged Product
data=V10m_forc; clear V10m_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

clear grdtype;

%% Downward LW:

grdtype='node';
var='lwr';

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/LW_forc_16.mat
data=LW_forc; clear LW_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Net Shortwave:

var='nswrs';

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/SW_forc_16.mat
data=SW_forc; clear SW_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Air - Pressure:

var='airpressure';

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/P_Surf_forc_16.mat
data=P_Surf_forc; clear P_Surf_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Surface air temperature

var='sat';

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/T_2m_air_forc_16.mat
data=T_2m_air_forc; clear T_2m_air_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Relative Humidity

var='rhumid';

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/RH_forc_16.mat
data=RH_forc; clear RH_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Precipitation 
% Store as Precipitation - Evaporation:

var = 'prate';

load D:/Studies/Atmospheric_Forcing/Merged_Product/2016/PE_forc_16.mat 
data=PE_forc; clear PE_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Evaporation
% Store zeros:

var = 'evap';
data=zeros(Mobj.nVerts,length(time));

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;