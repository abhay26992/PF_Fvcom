%% Script to obtain temperature @ the transect cell indices

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
%addpath(genpath('E:/Studies/Model_Out/Gamma_0.1/fvcom_4/Nest+Ice/PF_stdyr_runs/Output_years'));

%% 1. Input(s)

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Bathymetry_IceBridge/Input_Files/Natural_Sill/IBN_M.mat

my_case='qsd'; %or therm/dyn experiments
t_node=mk_ht_temp_ts_icediag(my_case); %load the temperature time-series 

%% 2a. Define the transect & get the indices 

% Note: Transects defined here should be consistent with those defined in
% mk_vtvars_ts_icediag.m to calculate the heat transport (u*T)

x1=[-2.846e+05 -2.717e+05];y1=[-9.313e+05 -9.197e+05]; %gyre-section

%x1=[-2.904e+05 -2.728e+05];y1=[-9.151e+05 -9.092e+05]; %fjord-mouth section
%x1=[-2.992e+05 -2.686e+05];y1=[-9.084e+05 -8.952e+05]; %fjord-mouth section - 1 

%% 2b. Calculate the indices

htvars.ind = fvcom_get_transect_ind(Mobj,x1,y1);

%% 3. Interpolate t_node to cells 

for i=1:size(t_node,3)
    t_cell(:,:,i)=node2cell(squeeze(t_node(:,:,i)),Mobj.tri);
end

htvars.temp=t_cell(htvars.ind,:,:);

htvars.temp_off=htvars.temp+2.5; %offset by 2.5 deg C
%% 4. Inspect

figure(1);clf
plot_field(Mobj,t_node(:,23,250));colormap jet;hold on;
scatter(Mobj.xc(htvars.ind),Mobj.yc(htvars.ind),50,squeeze(htvars.temp(:,23,250)));

%% 4. Save 

save htvars_qsd htvars  %or therm/dyn depending on my_case - gyre section

