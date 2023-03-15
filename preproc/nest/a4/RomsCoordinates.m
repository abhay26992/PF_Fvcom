function RCoor=RomsCoordinates(datasource,path)

%Reads ROMS coordinates from ....


%Read in from grid file
if strcmp(datasource,'nordic4km')
  %path='/home/olean/nordic4km/grid/';
  fname='ocean_grd.nc_nordic4km';
elseif strcmp(datasource,'norkyst800m')
  %path='/global/work/ovanov/NorKyst_800_Troms_2011/Results/';
  %fname='norkyst_800m_his_0200.nc';
  grdfile='/global/work/jsk/norkyst_800m_grid_full_WGS84.nc';
  fname='norkyst_800m_his.nc_2013122501-2013122600';
elseif strcmp(datasource,'arctic4_daily')
  grdfile='/global/work/apn/A4_modelruns/A4_nudging/linked_with_dates/a4_avg_2013-01-01.nc';
  fname='a4_avg_2013-01-01.nc';
end


ncid = netcdf.open([path fname],'NC_NOWRITE');


%lat and lon in rho points
id=netcdf.inqVarID(ncid,'lon_rho');
lon_rho=netcdf.getVar(ncid,id);lon_rho=double(lon_rho);
id=netcdf.inqVarID(ncid,'lat_rho');
lat_rho=netcdf.getVar(ncid,id);lat_rho=double(lat_rho);

%lat and lon i u,v points
  lon_u=0.5*(lon_rho(1:end-1,:)+lon_rho(2:end,:));
  lat_u=0.5*(lat_rho(1:end-1,:)+lat_rho(2:end,:));
  lon_v=0.5*(lon_rho(:,1:end-1)+lon_rho(:,2:end));
  lat_v=0.5*(lat_rho(:,1:end-1)+lat_rho(:,2:end));

%mask
id=netcdf.inqVarID(ncid,'mask_rho');
mask_rho=netcdf.getVar(ncid,id);mask_rho=double(mask_rho);
%id=netcdf.inqVarID(ncid,'mask_u');
%mask_u=netcdf.getVar(ncid,id);mask_u=double(mask_u);
%id=netcdf.inqVarID(ncid,'mask_v');
%mask_v=netcdf.getVar(ncid,id);mask_v=double(mask_v);
mask_u=ncread(grdfile,'mask_u');
mask_v=ncread(grdfile,'mask_v');
%pn and pm in rho points
id=netcdf.inqVarID(ncid,'pn');
pn=netcdf.getVar(ncid,id);pn=double(pn);
id=netcdf.inqVarID(ncid,'pm');
pm=netcdf.getVar(ncid,id);pm=double(pm);
dy=1./pn;
dx=1./pm;

%bottom topo
id=netcdf.inqVarID(ncid,'h');
h=netcdf.getVar(ncid,id);pm=double(pm);h=double(h);

%angle

id=netcdf.inqVarID(ncid,'angle');
roms_angle=netcdf.getVar(ncid,id);roms_angle=double(roms_angle);

if strcmp(datasource,'nordic4km')
    netcdf.close(ncid)
end


%%s-coordinate
if strcmp(datasource,'nordic4km')
  %Parameters from infile
  hc=10;%TCLINE
  theta_b=0.4;%THETA_B
  theta_s=3;%THETA_S
  N=32;
elseif strcmp(datasource,'arctic4_daily')
 id=netcdf.inqVarID(ncid,'hc');
 hc=netcdf.getVar(ncid,id);hc=double(hc);
 id=netcdf.inqVarID(ncid,'theta_b');
 theta_b=netcdf.getVar(ncid,id);theta_b=double(theta_b);
 id=netcdf.inqVarID(ncid,'theta_s');
 theta_s=netcdf.getVar(ncid,id);theta_s=double(theta_s);
 N=35;
 id=netcdf.inqVarID(ncid,'Vstretching');
 Vstretching=netcdf.getVar(ncid,id);
 id=netcdf.inqVarID(ncid,'Vtransform');
 Vtransform=netcdf.getVar(ncid,id);
elseif strcmp(datasource,'norkyst800m')
 id=netcdf.inqVarID(ncid,'hc');
 hc=netcdf.getVar(ncid,id);hc=double(hc);
 id=netcdf.inqVarID(ncid,'theta_b');
 theta_b=netcdf.getVar(ncid,id);theta_b=double(theta_b);
 id=netcdf.inqVarID(ncid,'theta_s');
 theta_s=netcdf.getVar(ncid,id);theta_s=double(theta_s);
 %dimid=netcdf.inqDimID(ncid,'N');
 dimid=netcdf.inqDimID(ncid,'s_rho');
 
 [tmp N]=netcdf.inqDim(ncid,dimid);
 id=netcdf.inqVarID(ncid,'Vstretching');
 Vstretching=netcdf.getVar(ncid,id);
 id=netcdf.inqVarID(ncid,'Vtransform');
 Vtransform=netcdf.getVar(ncid,id);
end


[Cs_rho s_rho]=roms_cs_calc(theta_s,theta_b,N,'r');
%Calculate z_rho and z_w, first initialize

 if strcmp(datasource,'nordic4km')
  nx=1024;ny=580;
  
elseif strcmp(datasource,'arctic4_daily')
  nx=1602;ny=1202;
 
elseif strcmp(datasource,'norkyst800m')
  dimid=netcdf.inqDimID(ncid,'xi_rho');
  [tmp nx]=netcdf.inqDim(ncid,dimid);
  dimid=netcdf.inqDimID(ncid,'eta_rho');
  [tmp ny]=netcdf.inqDim(ncid,dimid);
  
 end

 
 
 %netcdf.close(ncid)
%z_rho=NaN(nx,ny,N);
%z_w=NaN(nx,ny,N+1);

if strcmp(datasource,'norkyst800m')
zet=zeros(size(lon_rho));
elseif strcmp(datasource,'nordic4km')
load zeta_avg.mat
zet=zeros(size(zeta_avg));
elseif strcmp(datasource,'arctic4_daily') 
id=netcdf.inqVarID(ncid,'zeta');
zet=netcdf.getVar(ncid,id);
zet=squeeze(zet(:,:,1));
end



z_rho=set_depth(Vtransform, Vstretching, ...
                       theta_s, theta_b, hc, N, ...
		       1, h, zet, 0);
z_u=set_depth(Vtransform, Vstretching, ...
                       theta_s, theta_b, hc, N, ...
	               3, h, zet, 0);
z_v=set_depth(Vtransform, Vstretching, ...
                       theta_s, theta_b, hc, N, ...
	               4, h, zet, 0);
               
           
               
 RCoor.lon_rho=lon_rho;
 RCoor.lat_rho=lat_rho;
 %RCoor.lon_u=lon_u;
 %RCoor.lat_u=lat_u;
 %RCoor.lon_v=lon_v;
 %RCoor.lat_v=lat_v;
 RCoor.z_rho=z_rho;
 %RCoor.z_u=z_u;
 %RCoor.z_v=z_v;
 RCoor.h=h;
 RCoor.mask_rho=mask_rho;
 RCoor.mask_u=mask_u;
 RCoor.mask_v=mask_v;
 RCoor.roms_angle=roms_angle;
 
