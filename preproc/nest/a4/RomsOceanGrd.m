function ROGrd=RomsOceanGrd(RCoor,rind)

% RomsOceanGrd chooses ocean rho point based on the tailed grid which match
% the averaged u,v. It also generate index for finding the ocean 
% points 

% Input
%  RCoor: a structure which contains grid information for a chosen ROMS
%  it: index for tailoring the ROMS domain, which is determined manually


 use RCoor
 
 %%%%%crop norkyst domain to local domain
 
 use rind
 
 lon_rho=lon_rho(idr,jdr);
  lat_rho=lat_rho(idr,jdr); 
  z_rho=z_rho(idr,jdr,:);
  h=h(idr,jdr);
  mask_rho=mask_rho(idr,jdr);
  mask_u=mask_u(idu,jdu);
  mask_v=mask_v(idv,jdv);
  roms_angle=roms_angle(idr,jdr);
  
 %%%%%%%clip domain  to match averaged u and v%%%
 
  lon_rho=lon_rho(2:end-1,2:end-1);
  lat_rho=lat_rho(2:end-1,2:end-1); 
  z_rho=z_rho(2:end-1,2:end-1,:);
  h=h(2:end-1,2:end-1);
  mask_rho=mask_rho(2:end-1,2:end-1);
  mask_u=0.5.*(mask_u(1:end-1,:)+mask_u(2:end,:));
  mask_u=mask_u(:,2:end-1);
  mask_v=0.5.*(mask_v(:,1:end-1)+mask_v(:,2:end));
  mask_v=mask_v(2:end-1,:);
  roms_angle=roms_angle(2:end-1,2:end-1);
   
   
 
 
 ocn_r=find(mask_rho==1&mask_u==1&mask_v==1);
 
 
ell=referenceEllipsoid('wgs84','m');
a=ell.SemimajorAxis;e=ell.Eccentricity;
defaultellipsoid=[a,e];
[xr, yr, gam_abp, k_abp] = polarst_fwd(1, lat_rho, lon_rho, defaultellipsoid);

%[xr,yr]=deg2utm(lat_rho,lon_rho,'33 W'); % to match the fvcom u,v point : Edited -> May 8,2019

xr=xr(ocn_r);
yr=yr(ocn_r);

hr=h(ocn_r);

[nx,ny,nz]=size(z_rho);
zr=[];

for n=1:nz
    tmpz=squeeze(z_rho(:,:,n));
    zr(:,end+1)=tmpz(ocn_r);    
end


 ROGrd.xr=xr;
 ROGrd.yr=yr;
 ROGrd.hr=hr;
 ROGrd.zr=zr;
 ROGrd.ocn_r=ocn_r;
 ROGrd.angle=roms_angle;
 ROGrd.angler=roms_angle(ocn_r);
 
 
 
