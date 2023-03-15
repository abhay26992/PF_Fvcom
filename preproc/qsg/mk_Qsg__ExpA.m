%Build the river forcing file

%Written by Abhay Prakash (contact: abhay.prakash@natgeo.su.se)

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));
clear all;close all;clc;

% Experiment (a):
% In this experiment, we ignore the seasonality in the runoff data, and
% instead, take the summer mean value, and spread it uniformly across all
% the cells that make up the GL.

%% Input parameters 

load ./IBN_M.mat;
kbm1=size(Mobj.siglay,2);

%1. Time Component: Daily avg. fields
time=datenum(2016,01,01):datenum(2017,01,01); %Simulation Time
time=time-datenum(1858,11,17); %FVCOM time
nt=length(time);

%2. Space Component: Cell numbers of river mouth at the GL 
% Important Note: Find all cell numbers at the GL which would, together, 
% make up the subglacial discharge outlets. 

figure(1);clf
plot_field(Mobj,Mobj.zisf);colormap jet;caxis([0 500]);hold on;
%[xr,yr]=ginput; %x,y points along the GL
load ./xr;load ./yr;
ind_GL=find(inpolygon(Mobj.x,Mobj.y,xr,yr));hold on;
scatter(Mobj.x(ind_GL),Mobj.y(ind_GL),100,'k','o','filled');

%The tge.m script generates 'tge.mat' which contains 'isbce' - all cells 
%that fall on the boundary..and a lot more additional useful grid metrics -
%see tge.m script for details - for a change, it is nicely commented). 

%Mobj.tri = all nodes that make up a cell. Use this to find the cell #s 
%that contain our nodes (ismember). Then pass it through isbce to only 
%include those cells that are on the boundary.

load('./tge.mat')
bcell=find(isbce==1);%find all cells that are on a solid boundary (not obc)

glc_1=find(ismember(Mobj.tri(:,1),ind_GL));
glc_2=find(ismember(Mobj.tri(:,2),ind_GL));
glc_3=find(ismember(Mobj.tri(:,3),ind_GL));
glc=sort([glc_1;glc_2;glc_3]);

glc_bcell_ind=find(ismember(glc,bcell));
glc_bcell=unique(glc(glc_bcell_ind));
figure(2);clf
plot_field(Mobj,Mobj.zisf);colormap jet;caxis([0 500]);hold on;
scatter(Mobj.xc(glc_bcell),Mobj.yc(glc_bcell),100,'c','o','filled');

%Another approach: Find all cell indices (xc_ind,yc_ind) within a polygon
%region defined either using ginput or as 4 vertices (see subsetMobj.m 
%script). Note: Ensure that there are no lateral boundaries present here
%for isbce to work without any error. isbce can find cind in the interior
%of the mesh domain (0), on a solid boundary (1), open boundary (2), on 2
%solid boundaries (2)- therefore, if both lateral coastline and GL are
%present, isbce will include those indices as well.


%% Experiment A:

% Then,distribute the runoff over these points (can be done in any manner). 
% Here, without any information about the subglacial runoff routing, we
% distribute it evenly over the chosen points (Experiment A)

Q_sg=552; %Summer mean value in m^3/s, for now, taken from Mankoff_BMB

river.pnid=(glc_bcell)'; %Cell numbers of river mouth at the GL 
river.riverNames=num2str((1:numel(glc_bcell))'); %Corresponding river name

transport = repmat(Q_sg/numel(glc_bcell),[numel(glc_bcell),1]); %Equally distributed over all cells along the GL

cwdir=pwd;
cd C:\PhD\FVCOM\Matlab_Repository\Petermann_Bathy\Output\Gamma_0.1\fvcom_4\Nest+Atm+Ice\PF_stdyr_runs\Transport_Diagnostics;
fn='ocnice_2016_restart_0367.nc';
hc=ncread(fn,'h_center'); %depth @cell
cd(cwdir);
[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);

pval_hc=gsw_p_from_z(-hc(glc_bcell),latc(glc_bcell)); %pressure in dbar from h
temp=gsw_t_freezing(zeros(numel(glc_bcell),1),pval_hc); %pressure dependent freezing point of freshwater

salt=zeros(numel(glc_bcell),1); %salt =0

river.transport=repmat(transport,[1,nt]); % m^3/s
river.temp=repmat(temp,[1,nt]); % m^3/s
river.salt=repmat(salt,[1,nt]); % m^3/s
river.time=time;
river.vdist=ones(kbm1,1).*(1/kbm1); %vertical distribution - evenly in each layer [Note: sum of vdist must be 1]

%% Make the (separate) namelist and netcdf forcing file

use river

% Namelist:
np=length(pnid);

fname='riverNamelist.nml';
fid=fopen(fname,'w+');

riverFile='riverdata.nc';

formatspec = [repmat('%1.5f ',1,kbm1) '\n'];

for i=1:np
    
    fprintf(fid,' &NML_river\n');
    fprintf(fid,[' river_NAME = ' strcat(' ''',riverNames(i,:),''' ') '\n']);
    fprintf(fid,[' river_FILE = ' strcat(' ''',riverFile,''' ') '\n']);
    fprintf(fid,' river_GRID_LOCATION = ');
    fprintf(fid,'%6.0f\n',pnid(i));
    fprintf(fid,' RIVER_VERTICAL_DISTRIBUTION = ');
    fprintf(fid,formatspec,vdist);
    fprintf(fid,'/\n');
end
fclose(fid)


% Dumping river data (in river netcdf)


riverInfo1=[];
riverInfo2=[];
julian=1;

% river names (must be 80 character strings)
tmp  = '                                                                                ';
tmp=repmat(tmp,np,1);
tmpname=tmp';
for i=1:np
    fname =riverNames(i,:);
    tmpname(1:length(fname),i) = fname;
end

% write_FVCOM_river_Singlenode(riverFile,tmpname,time,transport,temp,salt,riverInfo1, riverInfo2, julian)

write_FVCOM_river_multinode(riverFile,tmpname,time,transport,temp,salt,riverInfo1, riverInfo2, julian)

