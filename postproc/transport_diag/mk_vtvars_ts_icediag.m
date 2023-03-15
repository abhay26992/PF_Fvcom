%% Get the u(x) and v(y) velocities at the transect indices

% Save all the intermediate variables & final product in a (vtvars) struct

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
%addpath(genpath('E:/Studies/Model_Out/Gamma_0.1/fvcom_4/Nest+Ice/PF_stdyr_runs'));

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Bathymetry_IceBridge/Input_Files/Natural_Sill/IBN_M.mat

%% 1a. Get cell indices along the transect joining (x1,y1) and (x2,y2)

%x1=[-2.846e+05 -2.717e+05];y1=[-9.313e+05 -9.197e+05]; %gyre-section
x1=[-2.825e+05 -2.663e+05];y1=[-9.506e+05 -9.394e+05]; %crack_calv-section

%x1=[-2.904e+05 -2.728e+05];y1=[-9.151e+05 -9.092e+05]; %fjord-mouth section
%x1=[-2.992e+05 -2.686e+05];y1=[-9.084e+05 -8.952e+05]; %fjord-mouth section - 1 

%% 1b. Calculate the indices

ind = fvcom_get_transect_ind(Mobj,x1,y1);
vtvars.ind=ind;

%figure(1);clf %view the section
%plot_field(Mobj,Mobj.zisf);colormap jet;caxis([0 400]);hold on;
%scatter(Mobj.xc(ind),Mobj.yc(ind),'r','o','filled');

%% 2. Get velocities at transect indices

my_case='qsd'; %or other qsg experiments
[vtvars.u,vtvars.v]=mk_transport_ts_icediag(vtvars.ind,my_case); %load the u & v time-series @ the transect indices

%% 3. Get anglec at transect indices

t_angle=transport_anglec(Mobj); %transport_anglec.m = changed utm2deg to polarstereo_inv in the anglec.m script
t_angle=t_angle(ind); % crop to transect
vtvars.anglec=repmat(t_angle,[1 size(Mobj.siglay,2) size(vtvars.u,3)]); % size = 180 x 23 x 366

%% 4. Vel_geo2xy using anglec

for i = 1:size(vtvars.u,3)
    [vtvars.ux(:,:,i),vtvars.vy(:,:,i)]=vel_geo2xy(vtvars.u(:,:,i),vtvars.v(:,:,i),vtvars.anglec(:,:,i));
end

%% 5. Save 

save vtvars_qsd vtvars  %or qsb/qsc etc. depending on my_case

