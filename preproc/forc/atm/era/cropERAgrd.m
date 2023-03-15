function erawdind=cropERAgrd(Mobj)

addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy'));
addpath(genpath('C:/PhD/FVCOM/Matlab_Repository/Stallo_Matlab'));

%Saves the original ERA grid in eragrd.mat and crops the wrf 
%grid around the fvcom grid given by Mobj. Find indices
%jdr, idr and stores them in the structure wdind.
lonmin=min(Mobj.lon)-0.1;lonmax=max(Mobj.lon)+0.1;
latmin=min(Mobj.lat)-0.1;latmax=max(Mobj.lat)+0.1;
%load roms grid
%fname='/global/work/hclear dj002/WRF3km/wrf_3km_norkyst800_2014.nc';
fname='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/ERA_Forcing/era_pmf_an_2010_19.nc';
lonr=ncread(fname,'longitude');
latr=ncread(fname,'latitude');

latr=latr';
[latera,lonera]=meshgrid(latr,lonr);

%Projection:
a = 6378137;
e = 0.08181919;
phi_c = 70;
lambda_0 = -45;
[xg, yg] = polarstereo_fwd(latera, lonera, a,e,phi_c,lambda_0);

save eragrd xg yg

%Find rows and colums in era grid where either lonr or latr 
%are within the fvcom limits
lontest=any(lonera>lonmin&lonera<lonmax&latera>latmin&latera<latmax);
erawdind.jdr=find(lontest);
lattest=any(latera>latmin&latera<latmax&lonera>lonmin&lonera<lonmax,2);
erawdind.idr=find(lattest);
figure(1);clf;
%pr=plot(latr(wdind.idr,wdind.jdr),lonr(wdind.idr,wdind.jdr),'.k');
pr=plot(lonera(erawdind.idr,erawdind.jdr),latera(erawdind.idr,erawdind.jdr),'.k');
hold on
p=plot(Mobj.lon,Mobj.lat,'.r');


save erawdind erawdind