
function rind=cropROMSgrd(Mobj)

addpath(genpath('/home/abhay/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('/home/abhay/Matlab_Repository/Stallo_Matlab'));

%Crops roms grid around the fvcom grid given by Mobj. Find indices
%of rho-point (ir,jr), v-point (iv,jv) and u-points (iu,ju), all
%stored in the structure rind.
lonmin=min(Mobj.lon);lonmax=max(Mobj.lon);
latmin=min(Mobj.lat);latmax=max(Mobj.lat);
%load roms grid
fname='/global/work/apn/A4_modelruns/A4_nudging/linked_with_dates/a4_avg_2013-01-01.nc';
lonr=ncread(fname,'lon_rho');
latr=ncread(fname,'lat_rho');
%Find rows and colums in roms grid where either lonr or latr 
%are within the fvcom limits
lontest=any(lonr>lonmin&lonr<lonmax&latr>latmin&latr<latmax);
rind.jdr=find(lontest);
lattest=any(latr>latmin&latr<latmax&lonr>lonmin&lonr<lonmax,2);
rind.idr=find(lattest);
%Find indices for u and v
rind.idu=rind.idr(2:end);rind.jdu=rind.jdr;
rind.idv=rind.idr;rind.jdv=rind.jdr(2:end);
save rind rind


clf,pr=plot(latr(rind.idr,rind.jdr),lonr(rind.idr,rind.jdr),'.k');
hold on
p=plot(Mobj.lat,Mobj.lon,'.r');
