function wdind=cropRACgrd(Mobj)

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));

%Saves the original RACMO grid in racgrd.mat and crops the rac 
%grid around the fvcom grid given by Mobj. Find indices
%jdr, idr and stores them in the structure wdind.
lonmin=min(Mobj.lon)-0.1;lonmax=max(Mobj.lon)+0.1;
latmin=min(Mobj.lat)-0.1;latmax=max(Mobj.lat)+0.1;
%load roms grid
%fname='/global/work/hclear dj002/WRF3km/wrf_3km_norkyst800_2014.nc';
fname='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/t2m.nc';
lonr=ncread(fname,'lon');
latr=ncread(fname,'lat');

%Projection:
a = 6378137;
e = 0.08181919;
phi_c = 70;
lambda_0 = -45;
[xg, yg] = polarstereo_fwd(latr, lonr, a,e,phi_c,lambda_0);

save racgrd xg yg

%Find rows and colums in roms grid where either lonr or latr 
%are within the fvcom limits
lontest=any(lonr>lonmin&lonr<lonmax&latr>latmin&latr<latmax);
wdind.jdr=find(lontest);
lattest=any(latr>latmin&latr<latmax&lonr>lonmin&lonr<lonmax,2);
wdind.idr=find(lattest);
clf;
%pr=plot(latr(wdind.idr,wdind.jdr),lonr(wdind.idr,wdind.jdr),'.k');
pr=plot(lonr(wdind.idr,wdind.jdr),latr(wdind.idr,wdind.jdr),'.k');
hold on
p=plot(Mobj.lon,Mobj.lat,'.r');


save wdind wdind