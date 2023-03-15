%% Merge and create the mean-monthly (2007-2009) climatology struct

% The 'standard-year' FVCOM mean-monthly climatology uses the A4 data from
% 2007-01-01 to 2009-01-01. Here, we merge the two structs month-wise which 
% were created using A4_clima.m and prepare the mean monthly climatology
% struct.

% Split in 3 parts. This is part (2) of (3). 

% See previous --> A4_clima.m (Part 1)
% See further --> A4_write_clima.m (Part 3)

% Written by Abhay Prakash (abhay.prakash@natgeo.su.se)

clear all;close all;clc;

%% Load climatologies

load ./A4_2007_2017/Boundary/2007-09/clima_07.mat;clima07=clima;clear clima;
load ./A4_2007_2017/Boundary/2007-09/clima_08.mat;clima08=clima;clear clima;

%% Concatenate month-wise and create mean (2007-2009) monthly climatology

for i=1:size(clima07.temp,3)
    temp(:,:,i)=mean(cat(3,clima07.temp(:,:,i),clima08.temp(:,:,i)),3);
    salt(:,:,i)=mean(cat(3,clima07.salt(:,:,i),clima08.salt(:,:,i)),3);
    u(:,:,i)=mean(cat(3,clima07.u(:,:,i),clima08.u(:,:,i)),3);
    v(:,:,i)=mean(cat(3,clima07.v(:,:,i),clima08.v(:,:,i)),3);
    ubar(:,i)=mean(cat(2,clima07.ubar(:,i),clima08.ubar(:,i)),2);
    vbar(:,i)=mean(cat(2,clima07.vbar(:,i),clima08.vbar(:,i)),2);
    zeta(:,i)=mean(cat(2,clima07.zeta(:,i),clima08.zeta(:,i)),2);
end

%% Save 

clima_0708.temp=temp;
clima_0708.salt=salt;
clima_0708.u=u;
clima_0708.v=v;
clima_0708.ubar=ubar;
clima_0708.vbar=vbar;
clima_0708.zeta=zeta;

save ./A4_2007_2017/Boundary/2007-09/clima_0708 clima_0708
