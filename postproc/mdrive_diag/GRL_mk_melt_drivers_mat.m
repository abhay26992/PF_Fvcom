%% Make melt + melt driver variables

% Can be extended to include additional experiments

clear all;close all;clc;
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));

load C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Bathymetry_IceBridge/Input_Files/Natural_Sill/IBN_M.mat

%% Output files

fn_std='D:/Abhay/Model_Out/PF_std_yr_2016/ocnice_2016_avg_0001.nc';
fn_qsa='D:/Abhay/Model_Out/PF_Q_sg.exp_a/exp_a_avg_0001.nc';
fn_qsb='D:/Abhay/Model_Out/PF_Q_sg.exp_b/exp_b_avg_0001.nc';
fn_qsc='D:/Abhay/Model_Out/PF_Q_sg.exp_c/exp_c_avg_0001.nc';
fn_qsd='D:/Abhay/Model_Out/PF_Q_sg.exp_d/exp_d_avg_0001.nc';

%% Load variables

%Standard Run
t_std=ncread(fn_std,'temp');t_std=squeeze(t_std(:,1,:));
s_std=ncread(fn_std,'salinity');s_std=squeeze(s_std(:,1,:));
u_std=ncread(fn_std,'u');u_std=squeeze(u_std(:,1,:));
v_std=ncread(fn_std,'v');v_std=squeeze(v_std(:,1,:));
melt_std=ncread(fn_std,'meltrate');

%Exp_A 
t_qsa=ncread(fn_qsa,'temp');t_qsa=squeeze(t_qsa(:,1,:));
s_qsa=ncread(fn_qsa,'salinity');s_qsa=squeeze(s_qsa(:,1,:));
u_qsa=ncread(fn_qsa,'u');u_qsa=squeeze(u_qsa(:,1,:));
v_qsa=ncread(fn_qsa,'v');v_qsa=squeeze(v_qsa(:,1,:));
melt_qsa=ncread(fn_qsa,'meltrate');

%Exp_B 
t_qsb=ncread(fn_qsb,'temp');t_qsb=squeeze(t_qsb(:,1,:));
s_qsb=ncread(fn_qsb,'salinity');s_qsb=squeeze(s_qsb(:,1,:));
u_qsb=ncread(fn_qsb,'u');u_qsb=squeeze(u_qsb(:,1,:));
v_qsb=ncread(fn_qsb,'v');v_qsb=squeeze(v_qsb(:,1,:));
melt_qsb=ncread(fn_qsb,'meltrate');

%Exp_C 
t_qsc=ncread(fn_qsc,'temp');t_qsc=squeeze(t_qsc(:,1,:));
s_qsc=ncread(fn_qsc,'salinity');s_qsc=squeeze(s_qsc(:,1,:));
u_qsc=ncread(fn_qsc,'u');u_qsc=squeeze(u_qsc(:,1,:));
v_qsc=ncread(fn_qsc,'v');v_qsc=squeeze(v_qsc(:,1,:));
melt_qsc=ncread(fn_qsc,'meltrate');

%Exp_D 
t_qsd=ncread(fn_qsd,'temp');t_qsd=squeeze(t_qsd(:,1,:));
s_qsd=ncread(fn_qsd,'salinity');s_qsd=squeeze(s_qsd(:,1,:));
u_qsd=ncread(fn_qsd,'u');u_qsd=squeeze(u_qsd(:,1,:));
v_qsd=ncread(fn_qsd,'v');v_qsd=squeeze(v_qsd(:,1,:));
melt_qsd=ncread(fn_qsd,'meltrate');


%% Call q_mk_fp

[Mobj.ntve,Mobj.nbve,Mobj.nbvt]=get_nbve(Mobj);

[Tdrive_std,Sdrive_std,ustar_std]=q_mk_fp(Mobj,s_std,t_std,u_std,v_std);
[Tdrive_qsa,Sdrive_qsa,ustar_qsa]=q_mk_fp(Mobj,s_qsa,t_qsa,u_qsa,v_qsa);
[Tdrive_qsb,Sdrive_qsb,ustar_qsb]=q_mk_fp(Mobj,s_qsb,t_qsb,u_qsb,v_qsb);
[Tdrive_qsc,Sdrive_qsc,ustar_qsc]=q_mk_fp(Mobj,s_qsc,t_qsc,u_qsc,v_qsc);
[Tdrive_qsd,Sdrive_qsd,ustar_qsd]=q_mk_fp(Mobj,s_qsd,t_qsd,u_qsd,v_qsd);

%% Save

save Tdrive_std Tdrive_std;save ustar_std ustar_std;save melt_std melt_std;
save Tdrive_qsa Tdrive_qsa;save ustar_qsa ustar_qsa;save melt_qsa melt_qsa; 
save Tdrive_qsb Tdrive_qsb;save ustar_qsb ustar_qsb;save melt_qsb melt_qsb;
save Tdrive_qsc Tdrive_qsc;save ustar_qsc ustar_qsc;save melt_qsc melt_qsc; 
save Tdrive_qsd Tdrive_qsd;save ustar_qsd ustar_qsd;save melt_qsd melt_qsd;


