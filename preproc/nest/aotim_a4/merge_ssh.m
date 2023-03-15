%% Merge SSH fields from AOTIM and ROMS for FVCOM

%Written by Abhay Prakash (abhay.prakash@natgeo.su.se)
%Revision History: Version 2 (14 June, 2019)

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
load ./M.mat
load ./ngrd.mat
%% Block 1: SSH from AOTIM and the first struct:

ncdisp('jul_el_obc.nc');
ssh_aotim=ncread('jul_el_obc.nc','elevation');
ssh_ao_time=ncread('jul_el_obc.nc','time'); %Modified Julian Day (MJD)

%Note: length(ssh_ao_time)=2161 hours from 1/06/06 to 29/08/06 (90 days) +
%30/08/06 00:00:00 (90*24 + 1 =2161 data points).


% Convert the time:
ssh_ao_time=double(ssh_ao_time); 
time_check_1=datestr(ssh_ao_time); %Currently in FVCOM time
ssh_ao_time=ssh_ao_time+678942; % Add the time offset from FVC to Matlab
time_check_2=datestr(ssh_ao_time); % Now in Matlab convention

%Add the time offset from AOTIM 2006 to ROMS 2016:
ssh_ao_time=ssh_ao_time+3653; %10 years (365.25 days / year)
time_check_3=datestr(ssh_ao_time); %Now starting at 01-June-2016 (to 30-Aug-2016)

%% Block 2: 

roms_zeta_time %Execute this script to get the ROMS SSH field (nclm_zeta) defined over the 91 day period
               %01/06/16 - 30/08/16 [same overlapping period as AOTIM], and the ROMS time field (nclm_time)
               %which are daily average (total=91) points.


%Part a: Get the lat,lon at FVCOM open boundary nodes

if Mobj.have_strings %Get the list of open boundary nodes.
    tmpObcNodes = Mobj.obc_nodes';
    oNodes = tmpObcNodes(tmpObcNodes ~= 0)';
    fvlon = Mobj.x(oNodes); fvlat = Mobj.y(oNodes); %Lat Lon indexed at oNodes  
    assert(length(fvlon) == length(fvlat), 'Inconsistent number of coordinates for the open boundary.')
    
ssh_aotim_romsa4=nan(size(ssh_aotim));    
    
    for p = 1:length(fvlon)
%Part b: Find the closest ROMS position and interpolate the SSH to the same
%sampling as the AOTIM SSH.
        fx = fvlon(p);fy = fvlat(p);
        %dist=sqrt((ngrd.xn(:) - fx).^2 + (ngrd.yn(:) - fy).^2); %Get the distance (between ngrd/ROMS lat,lon and (one node at a time) FVC obc node)
        %sort_dist=sort(dist); %Then sort (in ascending order/default); and lastly, get indices:
        [~, jj] = sort(sqrt((ngrd.xn(:) - fx).^2 + (ngrd.yn(:) - fy).^2));
        ssh = interp1(nclm_time, nclm_zeta(jj(1), :), ssh_ao_time'); %Interpolate
        ssh_aotim_romsa4(p, :) = ssh_aotim(p, :) + ssh; %Add the two ssh fields
    end
    %clear ssh ir ic cc fx fy p fvlon fvlat oNodes tmpObcNodes
end

% Plotting:
plot(ssh_ao_time,ssh_aotim(1,:),'r');hold on;
plot(ssh_ao_time,ssh_aotim_romsa4(1,:),'k');
datetick('x','yyyymm');legend('aotim','romsa4');
xlabel('Time');ylabel('Sea Surface Height [meters]');

%% Write out as netcdf file:
Mobj.surfaceElevation=ssh_aotim_romsa4;
ElevationFile='pf_n_to_jul_el_obc.nc';
MyTitle='PFjord 2D tides';

%Need to revert back to time_check_1 (FVCOM time) when writing the input
%file for the model run:
ssh_ao_time=ssh_ao_time-datenum(1858,11,17,0,0,0);
time_check_4=datestr(ssh_ao_time);

%Note: Our time is now in MJD [i.e. Julian Time - datenum(1858,11,17,0,0,0) = 157 years]
% or, days since 17-11-1858. We do need to change anything to obc_nodes.
% Therefore, we do not call set_elevtide.m and use write_FVCOM_elevtide
% directly.

write_FVCOM_elevtide(Mobj,ssh_ao_time,ElevationFile,MyTitle);