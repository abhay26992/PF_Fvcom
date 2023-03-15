%% Dumping the Atmospheric forcing files for FVCOM:

clear all;close all;clc;

%% Load common inputs to the write_FVCOM_forcing_var function 
load ./M 

load ./Atm_Forc_T_Vec.mat  % 3-hr time steps from RACMO
time = RACMO_time; clear RACMO_time; 
time=time-datenum(1858,11,17,0,0,0);

fileprefix='mergedforc';

%% U-V wind velocity:

%%%%%%%%%%%%%%%%%%%%%%%%%%%% U-10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grdtype='cell';
var='uwnd';

load ./U10m_forc.mat %The Merged Product
data=U10m_forc; clear U10m_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V-10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var='vwnd';

load ./V10m_forc.mat %The Merged Product
data=V10m_forc; clear V10m_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

clear grdtype;

%% Downward LW:

grdtype='node';
var='lwr';

load ./LW_forc.mat
data=LW_forc; clear LW_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Net Shortwave:

var='nswrs';

load ./SW_forc.mat
data=SW_forc; clear SW_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Air - Pressure:

var='airpressure';

load ./P_Surf_forc.mat
data=P_Surf_forc; clear P_Surf_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Surface air temperature

var='sat';

load ./T_2m_air_forc.mat
data=T_2m_air_forc; clear T_2m_air_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Relative Humidity

var='rhumid';

load ./RH_forc.mat
data=RH_forc; clear RH_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Precipitation 
% Store as Precipitation - Evaporation:

var = 'prate';

load ./PE_forc.mat 
data=PE_forc; clear PE_forc;

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Evaporation
% Store zeros:

var = 'evap';
data=zeros(Mobj.nVerts,length(time));

write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Cloud - Cover (for A4-CICE)

var= 'cloudcover';

load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/TCC_ERA_2016.mat
data = TCC_ERA; era_time=time_era; clear TCC_ERA; clear time_era;

load ./PF_Multiyear_Run/Atm_Forc_T_Vec_2016.mat
hs_data=zeros(length(data),length(RACMO_time));

for p=1:length(data)
    hs_data(p,:)=interp1(era_time',data(p,:),RACMO_time');
end
clear data; clear era_time;


time=RACMO_time; clear RACMO_time;
time=time-datenum(1858,11,17,0,0,0);

data=hs_data; clear hs_data;


write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;

%% Specific Humidity (for A4-CICE)

var = 'spq';

load D:/Studies/Atmospheric_Forcing/ERA/ERA_Forcing/2016/SH_ERA_2016.mat
data=SH_ERA; era_time=time_era; clear SH_ERA; clear time_era;

load ./PF_Multiyear_Run/Atm_Forc_T_Vec_2016.mat
hs_data=zeros(length(data),length(RACMO_time));

for p=1:length(data)
    hs_data(p,:)=interp1(era_time',data(p,:),RACMO_time');
end
clear data; clear era_time;


time=RACMO_time; clear RACMO_time;
time=time-datenum(1858,11,17,0,0,0);

data=hs_data; clear hs_data;


write_FVCOM_forcing_var(Mobj, fileprefix, var,data, time);
clear var;clear data;