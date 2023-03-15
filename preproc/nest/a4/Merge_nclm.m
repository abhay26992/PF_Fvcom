%% Merge all the ROMS data structs to be written out in the nesting file

% Written by Abhay Prakash (abhay.prakash@natgeo.su.se)

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
% 
clear all;close all;clc;
%% 

%year=['2014'];month=['07';'08';'09';'10';'11';'12'];
year=['2014'];month=['01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'];
gen_name='pf_o_1_nest_';

for i=1:size(year,1)
    for j=1:size(month,1)
        filename=[gen_name, year(i,:),'_',month(j,:),'_all_days.mat'];
        yr=year(i,:);mnt=month(j,:);
        nclm=importdata(filename);
        eval(['nclm' yr(i,:) mnt(i,:) '= nclm']);
        clear nclm;
    end
end

% Manipulate to add Jan - 01 - YYYY - 00:00:00 (from the next year)

load pf_o_1_nest_2015_01_all_days.mat; nclm201501=nclm; clear nclm


%% Merge the nclm structs into one struct:

%Concatenate across dimension=time

%nclm=CatStructFields(nclm201407,nclm201408,nclm201409,nclm201410,nclm201411,nclm201412,nclm201501);
nclm=CatStructFields(nclm201401,nclm201402,nclm201403,nclm201404,nclm201405,nclm201406,nclm201407,nclm201408,nclm201409,nclm201410,nclm201411,nclm201412,nclm201501);

%Remove last 30 entries that correspond to January 2 - January 31

%t_stamp=nclm.time;check_t_stamp=datestr(t_stamp);
%new_t_stamp=t_stamp(1:end-30);check_new_t_stamp=datestr(new_t_stamp);

nclm.u=nclm.u(:,:,1:end-30);nclm.v=nclm.v(:,:,1:end-30);
nclm.salt=nclm.salt(:,:,1:end-30);nclm.temp=nclm.temp(:,:,1:end-30);
nclm.ubar=nclm.ubar(:,1:end-30);nclm.vbar=nclm.vbar(:,1:end-30);
nclm.zeta=nclm.zeta(:,1:end-30);nclm.time=nclm.time(:,1:end-30);

%Save:
save nclm_14 nclm



%% Write out to the netcdf nesting file:

clear all;close all;clc;

load ./ngrd; load ./nclm_14; time=nclm.time; nest_type=3; geo2xy=0;
fileprefix='pf_jan_dec_14';

write_FVCOM_nest(ngrd, fileprefix, nclm,time,nest_type,geo2xy)
%% Check nesting:
fn='pf_jan_dec_14.nc';
load ./M.mat;

%Plot the model bathymetry:
figure(1);clf
plot_field(Mobj,Mobj.h,'coordinate','spherical');colormap jet;

nestlat=ncread(fn,'lat');
nestlon=ncread(fn,'lon');
nesth=ncread(fn,'h');

%

hold on;%Overlay the bathymetry from the nested setup
scatter(nestlon,nestlat,100,nesth)
