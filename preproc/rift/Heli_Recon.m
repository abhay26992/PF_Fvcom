%% Helicopter Recon

% Check the coordinates of the Petermann Crack from the Helicopter survey
% and how it compares with your Sentinel-2A 10-m digitised crack.

%% The Sentinel Crack

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;

Crack=csvread('Crack_XY.csv',1,0);
x=Crack(:,1);y=Crack(:,2);

UTM_Zone='21 N';
UTM_Zone=repmat(UTM_Zone,length(x),1);
[Lat,Lon]=utm2deg(x,y,UTM_Zone);

%% Helicopter Crack
Heli_Crack=csvread('PetermannCrack_20190904_helitrack.csv',18,0);
Heli_Lon=Heli_Crack(:,1);Heli_Lat=Heli_Crack(:,2);

load ./M.mat

% Visualise the crack over 'zisf'
figure(1);clf
plot_field(Mobj,Mobj.zisf,'coordinate','spherical');colormap jet;
caxis([0 400]);
hold on;
scatter(Heli_Lon,Heli_Lat,20,'k','o','filled'); % Helicopter 
scatter(-61.373,80.84,'r','o','filled'); %The spot location from Martin
scatter(Lon,Lat,30,'y','*'); % Sentinel


