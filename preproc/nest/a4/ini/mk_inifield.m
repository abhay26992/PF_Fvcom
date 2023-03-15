%mk_inifield makes initial field for making a restart file for the FVCOM
%run

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
romsfile='./a4_avg_2016-06-01.nc';

ncid = netcdf.open(romsfile,'NC_NOWRITE');

tid=1; % 00:00:00 of the next day
%zeta
zeta=ncread(romsfile,'zeta',[idr(1) jdr(1) tid],[irl jrl 1]);
zeta=zeta(2:end-1,2:end-1,:);
%salt
salt=ncread(romsfile,'salt',[idr(1) jdr(1) 1 tid],[irl jrl Inf 1]);
salt=salt(2:end-1,2:end-1,:,:);
%temp
temp=ncread(romsfile,'temp',[idr(1) jdr(1) 1 tid],[irl jrl Inf 1]);
temp=temp(2:end-1,2:end-1,:,:);
%u
u=ncread(romsfile,'u',[idu(1) jdu(1) 1 tid],[iul jul Inf 1]);
u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:));
u=u(:,2:end-1,:,:);
%ubar
ubar=ncread(romsfile,'ubar',[idu(1) jdu(1) tid],[iul jul 1]);
ubar=0.5*(ubar(1:end-1,:,:)+ubar(2:end,:,:));
ubar=ubar(:,2:end-1,:);
%v
v=ncread(romsfile,'v',[idv(1) jdv(1) 1 tid],[ivl jvl Inf 1]);
v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));
v=v(2:end-1,:,:,:);
%vbar
vbar=ncread(romsfile,'vbar',[idv(1) jdv(1) tid],[ivl jvl 1]);
vbar=0.5*(vbar(:,1:end-1,:)+vbar(:,2:end,:));
vbar=vbar(2:end-1,:,:);

netcdf.close(ncid)


%transform u, v to true east and noth

[ubargeo,vbargeo]=vel_xy2geo(ubar,vbar,angle);



id=1;
%size(lon)
%size(lat)
%size(zeta)
figure,pcolor(lon,lat,zeta);shading interp
colorbar
caxis([-0.5 0.5])
hold on
quiver(lon(1:id:end,1:id:end),lat(1:id:end,1:id:end),ubargeo(1:id:end,1:id:end),vbargeo(1:id:end,1:id:end),10,'k');

[nx,ny,nz]=size(u);

ugeo=NaN(nx,ny,nz);
vgeo=NaN(nx,ny,nz);



for i=1:nz
    
    
    [ugeo(:,:,i),vgeo(:,:,i)]=vel_xy2geo(squeeze(u(:,:,i)),squeeze(v(:,:,i)),angle);
    
    
end


%interpolate zeta, ubar and var%



%zeta

Z0=zeta(ocn_r);

ini.zeta=sum(coef_r.*Z0(ind_r),2);



%ubar

UB0=ubargeo(ocn_r);

UBAR=sum(coef_r2c.*UB0(ind_r2c),2);

%var


VB0=vbargeo(ocn_r);

VBAR=sum(coef_r2c.*VB0(ind_r2c),2);


% convert to fvcom cartisian coordiante
fvangle=anglec(Mobj);
[UBAR1,VBAR1]=vel_geo2xy(UBAR,VBAR,fvangle);

%%
ini.ubar=UBAR;
ini.vbar=VBAR;

clear zeta ubar vbar

%%



%%%%%%%%%%%%%%%%%%%interpolate temperature,salinity and 3D velocities %%%%%%%%%%%%%%%%%%%%%%%

[nx,ny,nz]=size(salt);


% horizontal interpolation
disp('horizontal interpolation')
for i=1:nz
    
    %salt
    
    S0=squeeze(salt(:,:,i));
    
    S0=S0(ocn_r);
    
    SALT(:,i)=sum(coef_r.*S0(ind_r),2);
    
    %temp
    
    T0=squeeze(temp(:,:,i));
    
    T0=T0(ocn_r);
    
    TEMP(:,i)=sum(coef_r.*T0(ind_r),2);
    
    %u
    
    U0=squeeze(ugeo(:,:,i));
    
    U0=U0(ocn_r);
    
    U(:,i)=sum(coef_r2c.*U0(ind_r2c),2);
    
    %v
    
    V0=squeeze(vgeo(:,:,i));
    
    V0=V0(ocn_r);
    
    V(:,i)=sum(coef_r2c.*V0(ind_r2c),2);
    
    
    
    
end

clear salt temp u v ugeo vgeo

%%
%%%vertical interpolation
disp(['vertical interpolation'] )


salt=NaN(nn,kb-1);
temp=NaN(nn,kb-1);
u=NaN(nc,kb-1);
v=NaN(nc,kb-1);



for i=1:nn
    
    coef=squeeze(zcoef_n(i,:,:));
    ind=squeeze(zind_n(i,:,:));
    
    
    S1=squeeze(SALT(i,:));
    T1=squeeze(TEMP(i,:));
    
    salt(i,:)=sum(coef.*S1(ind),2);
    temp(i,:)=sum(coef.*T1(ind),2);
    
    clear coef ind
    
end

for i=1:nc
    
    coef=squeeze(zcoef_c(i,:,:));
    ind=squeeze(zind_c(i,:,:));
    U1=squeeze(U(i,:));
    V1=squeeze(V(i,:));
    u(i,:)=sum(coef.*U1(ind),2);
    v(i,:)=sum(coef.*V1(ind),2);
    
    
end

%%%% convert to the fvcom coordinate system

%[u,v]=vel_geo2xy(u,v,repmat(fvangle,1,kb-1));


clear U V



ini.salt=salt;
ini.temp=temp;
ini.u=u;
ini.v=v;

ini0=ini;

%%%%%%%%%%%%%%%%%%%%%%%smoothing zeta, salinity, temperature

ini.zeta=smoothfield(ini.zeta,Mobj,1,6);

%%%salt and temp

for i=1:kb-1
    ini.temp(:,i)=smoothfield(ini.temp(:,i),Mobj,1,6);
    ini.salt(:,i)=smoothfield(ini.salt(:,i),Mobj,1,6);
end

%%need to change the velocity back to cartisian with geo2xy

%save ini2014071001 ini ini0
save ini20160601 ini ini0

