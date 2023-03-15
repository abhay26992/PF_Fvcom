

clear all
close all
fabm=false;
zerovel=true;
saltcorr=false;
if saltcorr
  saltmin=10;
end

%load ini2014071001
load ini20160601
use ini0

if zerovel
 u=zeros(size(u));
 v=zeros(size(v));
 ubar=zeros(size(ubar));
 vbar=zeros(size(vbar));
end

if saltcorr
 i=find(salt<saltmin);
 salt(i)=saltmin;
end

load ./M
use Mobj
 

%Write data to restartup netcdf file.
casename=get_cstr;

%fname=['input/' casename '_restart_.nc'];
%fname='./input/Restart_2014071001.nc';
fname='pf_ini_rt_restart_0001.nc'
ncid=netcdf.open(fname,'nc_write');

tempid=netcdf.inqVarID(ncid,'temp');
saltid=netcdf.inqVarID(ncid,'salinity');
zetaid=netcdf.inqVarID(ncid,'zeta');
uid=netcdf.inqVarID(ncid,'u');
vid=netcdf.inqVarID(ncid,'v');
ubarid=netcdf.inqVarID(ncid,'ua');
vbarid=netcdf.inqVarID(ncid,'va');
%wid=netcdf.inqVarID(ncid,'ww');
etid=netcdf.inqVarID(ncid,'et');
rhoid=netcdf.inqVarID(ncid,'rho1');

if fabm
 load fabm_ini
 for ifabm=1:length(fini)
  eval(['fabm' num2str(ifabm) 'id=netcdf.inqVarID(ncid,fini(' num2str(ifabm) ').varname);'])
 end
end

%Initialize
%zeta=zeros(size(zeta));

netcdf.putVar(ncid,tempid,temp);
netcdf.putVar(ncid,saltid,salt);

netcdf.putVar(ncid,zetaid,zeta);
netcdf.putVar(ncid,uid,u);
netcdf.putVar(ncid,vid,v);
netcdf.putVar(ncid,ubarid,ubar);
netcdf.putVar(ncid,vbarid,vbar);
%netcdf.putVar(ncid,wid,w);
netcdf.putVar(ncid,etid,zeta);
rho=rho_F(salt,temp,-Mobj.siglayz);
netcdf.putVar(ncid,rhoid,rho);

if fabm
for ifabm=1:length(fini)
 eval(['netcdf.putVar(ncid,fabm' num2str(ifabm) 'id,fini(' num2str(ifabm) ').data);'])
end
end
%Initialize all other variables that have time as dimesion. FOr now we try to set them to zero.

netcdf.close(ncid);


