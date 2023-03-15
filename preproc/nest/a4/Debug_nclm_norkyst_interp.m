%% Debug nclm_norkyst_interp.m

load ./M.mat;

% Horizontal Interp
figure(1);clf;
plot_field(Mobj,Mobj.h);colormap jet;
hold on;scatter(xr,yr,100,T0); %T0,S0
caxis([-2 8]);

figure(2);clf;
plot_field(Mobj,Mobj.h);hold on;
scatter(lonn,latn,100,TT(:,1)); %SS,TT
colormap jet;caxis([-2 8]);colorbar;

A_=find(isnan(T0)); %T0,S0
T0(nind_r);ncoef_r.*ans;sum(ans,2);find(~isnan(ans)); %T0,S0


% Vertical Interp
load ./M.mat;
figure(3);clf;plot_field(Mobj,Mobj.h);hold on;scatter(xn,yn,100,temp(:,1)); %temp,salt
colormap jet;caxis([-2 8]);


%FVCOM on top of ROMS
figure(6);pcolor(lon_rho,lat_rho,temp_r(:,:,1));shading flat;hold on;scatter(lonn,latn,100,temp(:,23)) %ROMS bottom (1) on FVC bottom (23)
figure(7);pcolor(lon_rho,lat_rho,temp_r(:,:,35));shading flat;hold on;scatter(lonn,latn,100,temp(:,1)) %ROMS surface (35) on FVC surface (1)


%ROMS Horizontal Interp product (UU) over romsfile u
figure(7);pcolor(lon_rho,lat_rho,u(:,:,35));shading flat;hold on;scatter(lonc,latc,100,UU(:,35)) %ROMS surface = 35
figure(8);pcolor(lon_rho,lat_rho,u(:,:,1));shading flat;hold on;scatter(lonc,latc,100,UU(:,1)); %ROMS bottom =1

romspath='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Roms_a4_daily_nest/FV_R_PFjord/';
ncname='a4_avg_';year=['2016'];month=['06';'07';'08'];
load('pf_o_1_nest_2016_06_all_days.mat');load('ROGrd.mat')
nclm_norkyst_interp(romspath,ncname,year,month)


%% Nesting nc file -> Check NaNs

fn='pf_new_topo_nest.nc';
ncdisp(fn)
nv=ncread(fn,'nv');hyw=ncread(fn,'hyw');sal=ncread(fn,'salinity');
temp=ncread(fn,'temp');v=ncread(fn,'v');
u=ncread(fn,'u');ua=ncread(fn,'ua');va=ncread(fn,'va');
levc=ncread(fn,'siglev_center');layc=ncread(fn,'siglay_center');
lev=ncread(fn,'siglev');lay=ncread(fn,'siglay');
hcen=ncread(fn,'h_center');h=ncread(fn,'h');zeta=ncread(fn,'zeta');
lat=ncread(fn,'lat');lon=ncread(fn,'lon');x=ncread(fn,'x');y=ncread(fn,'y');
latc=ncread(fn,'latc');lonc=ncread(fn,'lonc');xc=ncread(fn,'xc');yc=ncread(fn,'yc');

find(isnan(zeta));
find(isnan(yc));
find(isnan(y));
find(isnan(xx));
find(isnan(xc));
find(isnan(x));
find(isnan(va));
find(isnan(v));
find(isnan(ua));
find(isnan(u));
find(isnan(temp));
find(isnan(sal));
find(isnan(nv));
find(isnan(lonc));
find(isnan(lon));
find(isnan(levc));
find(isnan(lev));
find(isnan(layc));
find(isnan(lay));
find(isnan(latc));
find(isnan(lat));
find(isnan(hyw));
find(isnan(hcen));
find(isnan(h));

% Plot:
load ./M.mat; x=ncread(fn,'x');y=ncread(fn,'y');h=ncread(fn,'h');
figure(3);clf
plot_field(Mobj,Mobj.h);colormap jet;hold on;scatter(x,y,100,h);

%% Roms File Output: Check for high surface temperature values

yy=ncread('a4_avg_2016-06-02.nc','temp');
xx=ncread('a4_avg_2016-06-20.nc','temp');
aa=ncread('a4_avg_2016-07-20.nc','temp');


zz=xx-yy; % day 20 - day 2
bb=aa-yy; % day 51 - day 2

figure(4);pcolorjw(zz(:,:,1));colorbar;colormap jet;
figure(5);pcolorjw(bb(:,:,1));colorbar;colormap jet;


% Should get progressively warmer if the results are correct
figure(6);pcolorjw(yy(:,:,1));colorbar;colormap jet;caxis([0 4]);
figure(7);pcolorjw(xx(:,:,1));colorbar;colormap jet;caxis([0 4]); 
figure(8);pcolorjw(aa(:,:,1));colorbar;colormap jet;caxis([0 4]);

%% Nesting File contents:

load ./M.mat;
fn='pf_new_topo_nest.nc';
temp=ncread(fn,'temp');
lat=ncread(fn,'lat');lon=ncread(fn,'lon');x=ncread(fn,'x');y=ncread(fn,'y');

figure(9);clf
%plot_field(Mobj,Mobj.h);hold on;
scatter(lon,lat,100,temp(:,1,20));colorbar;colormap jet;caxis([-2 7.5]);

% Restart file contents:
fn1='pf_ini_rt_restart_0001.nc';
r_temp=ncread(fn1,'temp');
rlat=ncread(fn1,'lat');rlon=ncread(fn1,'lon');rx=ncread(fn1,'x');ry=ncread(fn1,'y');

figure(10);clf
scatter(rlon,rlat,100,r_temp(:,1));colorbar;colormap jet;caxis([-2 4]);
