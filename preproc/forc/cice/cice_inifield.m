%mk_inifield makes initial field for making a restart file for the FVCOM
%run

%addpath(genpath('/home/abhay/Matlab_Repository/Petermann_Bathy'));
%addpath(genpath('/home/abhay/Matlab_Repository/Stallo_Matlab'));


clear all

close all

load ./M

xn=Mobj.x;
yn=Mobj.y;
xc=Mobj.xc;
yc=Mobj.yc;
nn=length(xn);
nc=length(xc);
kb=get_kb;

load rind
use rind
irl=length(idr);
jrl=length(jdr);

iul=length(idu);
jul=length(jdu);

ivl=length(idv);
jvl=length(jdv);

load ROGrd

use ROGrd  % cid fvangle h latc latn lonc lonn nid nv R cc xn yc yn

load RCoor

lat=RCoor.lat_rho(idr,jdr);
lon=RCoor.lon_rho(idr,jdr);
lat=lat(2:end-1,2:end-1);
lon=lon(2:end-1,2:end-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load nearest4_r  % coef_r ind_r
load nearest4_r2c
load Vinterp     % ZN ZC zind_c zind_n
use Vinterp


%%%%%%%%%%%read fields from roms output%%%%%%%%%%%%%%%%%%

% romsfile= '/global/work/apn/NK800/2014/norkyst_800m_his.nc_2014070901-2014071000';
%romsfile= '/global/work/apn/NK800/2014/norkyst_800m_his.nc_2014071001-2014071100';
%romsfile='/global/work/apn/A4_modelruns/A4_nudging/cice/iceh.2015-01-01.nc';
romsfile='C:\PhD\FVCOM\Matlab_Repository\Petermann_Bathy\Setup_Nesting\Roms_a4_daily_nest\FV_R_PFjord\iceh.2014-07-01.nc';

ncid = netcdf.open(romsfile,'NC_NOWRITE');

tid=1; % 00:00:00 of the next day

%aicen
aicen_d=ncread(romsfile,'aicen_d',[idr(1) jdr(1) 1 tid],[irl jrl Inf 1]);
aicen_d=aicen_d(2:end-1,2:end-1,:,:);

%vicen
vicen_d=ncread(romsfile,'vicen_d',[idr(1) jdr(1) 1 tid],[irl jrl Inf 1]);
vicen_d=vicen_d(2:end-1,2:end-1,:,:);

%vsnon
vsnon_d=ncread(romsfile,'vsnon_d',[idr(1) jdr(1) 1 tid],[irl jrl Inf 1]);
vsnon_d=vsnon_d(2:end-1,2:end-1,:,:);

%sice
sice_d=ncread(romsfile,'sice_d',[idr(1) jdr(1) tid],[irl jrl 1]);
sice_d=sice_d(2:end-1,2:end-1,:);

%tsfc
Tsfc_d=ncread(romsfile,'Tsfc_d',[idr(1) jdr(1) tid],[irl jrl 1]);
Tsfc_d=Tsfc_d(2:end-1,2:end-1,:);


%uvel
uvel_d=ncread(romsfile,'uvel_d',[idu(1) jdu(1) tid],[iul jul 1]);
uvel_d=0.5*(uvel_d(1:end-1,:,:)+uvel_d(2:end,:,:));
uvel_d=uvel_d(:,2:end-1,:);

%vvel
vvel_d=ncread(romsfile,'vvel_d',[idv(1) jdv(1) tid],[ivl jvl 1]);
vvel_d=0.5*(vvel_d(:,1:end-1,:)+vvel_d(:,2:end,:));
vvel_d=vvel_d(2:end-1,:,:);

netcdf.close(ncid)


%transform u, v to true east and noth

[uvelgeo,vvelgeo]=vel_xy2geo(uvel_d,vvel_d,angle);



id=1;
%size(lon)
%size(lat)
%size(zeta)
figure,pcolor(lon,lat,Tsfc_d);shading interp
colorbar
caxis([-0.5 0.5])
hold on
quiver(lon(1:id:end,1:id:end),lat(1:id:end,1:id:end),uvelgeo(1:id:end,1:id:end),vvelgeo(1:id:end,1:id:end),10,'k');

%[nx,ny,nz]=size(u);

%ugeo=NaN(nx,ny,nz);
%vgeo=NaN(nx,ny,nz);



%for i=1:nz
    
    
%    [ugeo(:,:,i),vgeo(:,:,i)]=vel_xy2geo(squeeze(u(:,:,i)),squeeze(v(:,:,i)),angle);
    
    
%end


%interpolate zeta, ubar and var%



%sice

SI=sice_d(ocn_r);
ini.sice_d=sum(coef_r.*SI(ind_r),2);


%Tsfc
T0=Tsfc_d(ocn_r);
ini.Tsfc_d=sum(coef_r.*T0(ind_r),2);


%uvel

UB0=uvelgeo(ocn_r);
UVEL=sum(coef_r2c.*UB0(ind_r2c),2);


%vvel


VB0=vvelgeo(ocn_r);
VVEL=sum(coef_r2c.*VB0(ind_r2c),2);


% convert to fvcom cartisian coordiante
fvangle=anglec(Mobj);
[UVEL1,VVEL1]=vel_geo2xy(UVEL,VVEL,fvangle);

%%
ini.uvel_d=UVEL;
ini.vvel_d=VVEL;

clear sice_d Tsfc_d uvel_d vvel_d

%%



%%%%%%%%%%%%%%%%%%%interpolate temperature,salinity and 3D velocities %%%%%%%%%%%%%%%%%%%%%%%

[nx,ny,nz]=size(aicen_d);


% horizontal interpolation
disp('horizontal interpolation')
for i=1:nz
    
    %aicen    
    AI=squeeze(aicen_d(:,:,i));
    
    AI=AI(ocn_r);
    
    AICEN_D(:,i)=sum(coef_r.*AI(ind_r),2);
    

    %vicen    
    VI=squeeze(vicen_d(:,:,i));
    
    VI=VI(ocn_r);
    
    VICEN_D(:,i)=sum(coef_r.*VI(ind_r),2);
    

    %vsnon 
    VS=squeeze(vsnon_d(:,:,i));
    
    VS=VS(ocn_r);
    
    VSNON_D(:,i)=sum(coef_r.*VS(ind_r),2);
    
    
    
end

clear aicen_d vicen_d vsnon_d

%%
%%%vertical interpolation
%disp(['vertical interpolation'] )


%salt=NaN(nn,kb-1);
%temp=NaN(nn,kb-1);
%u=NaN(nc,kb-1);
%v=NaN(nc,kb-1);



%for i=1:nn
    
%    coef=squeeze(zcoef_n(i,:,:));
%    ind=squeeze(zind_n(i,:,:));
    
    
%    S1=squeeze(SALT(i,:));
%    T1=squeeze(TEMP(i,:));
    
%    salt(i,:)=sum(coef.*S1(ind),2);
%    temp(i,:)=sum(coef.*T1(ind),2);
    
%    clear coef ind
    
%end

%for i=1:nc
    
%    coef=squeeze(zcoef_c(i,:,:));
%    ind=squeeze(zind_c(i,:,:));
%    U1=squeeze(U(i,:));
%    V1=squeeze(V(i,:));
%    u(i,:)=sum(coef.*U1(ind),2);
%    v(i,:)=sum(coef.*V1(ind),2);
    
    
%end

%%%% convert to the fvcom coordinate system

%[u,v]=vel_geo2xy(u,v,repmat(fvangle,1,kb-1));


%clear U V



ini.aicen_d=AICEN_D;
ini.vicen_d=VICEN_D;
ini.vsnon_d=VSNON_D;

ini0=ini;

%%%%%%%%%%%%%%%%%%%%%%%smoothing zeta, salinity, temperature

%ini.zeta=smoothfield(ini.zeta,Mobj,1,6);

%%%salt and temp

%for i=1:kb-1
%    ini.temp(:,i)=smoothfield(ini.temp(:,i),Mobj,1,6);
%    ini.salt(:,i)=smoothfield(ini.salt(:,i),Mobj,1,6);
%end

%%need to change the velocity back to cartisian with geo2xy

%save ini2014071001 ini ini0
save /home/abhay/Matlab_Repository/Petermann_Bathy/Setup_Nesting/fv_tools/FV_CICE_PFjord/Init_A4CICE/Init_A4CICE_2010_2018/ini20150101 ini ini0

