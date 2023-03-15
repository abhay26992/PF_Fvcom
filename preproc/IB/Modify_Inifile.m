% Modify the geometry (h and zisf) in the A4 initial condition file:

% Script written by Abhay Prakash (abhay.prakash@natgeo.su.se)
%% 

clear all;close all;clc;

%% The current inifile:

fn='ocnice_2015_restart_0001.nc';
ncdisp(fn)
h=ncread(fn,'h');

figure(1);clf
load ./M.mat
plot_field(Mobj,h);colormap jet;caxis([0 900]);

%% The new geometry:

clear Mobj; load ./Input_Files/Synthetic_Sill/IB_Exaggerated/IBE_M.mat;
% Check to confirm
figure(2);clf
plot_field(Mobj,Mobj.h);colormap jet;caxis([0 900]);

%% Modify the inifile:

ncid=netcdf.open(fn,'nc_write');

h_id=netcdf.inqVarID(ncid,'h');
netcdf.putVar(ncid,h_id,Mobj.h);

netcdf.close(ncid);

% Inspect:
figure(3);clf
changed_h=ncread(fn,'h');
plot_field(Mobj,changed_h);colormap jet;caxis([0 900]);

figure(4);clf
zisf=ncread(fn,'zisf');changed_wct=changed_h-zisf;
plot_field(Mobj,changed_wct);colormap jet;caxis([0 600]);