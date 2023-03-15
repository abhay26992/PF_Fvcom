clear all;close all;clc;
year=2016;

load ./wdind  % only uses the racmo data within the model domain
use wdind
irl=length(idr);
jrl=length(jdr);


% interpolating index and coefficients
load ./N4RAC
use N4RAC

nn=length(nindex);
nc=length(cindex);

%% Read in RACMO data, pay attention to the signs and units%%%%%%%%%%%

%% U-component (in m/s)

folder='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/';
file='u10m.nc';


% Change time for RACMO (Days since "1950-01-01")
time=ncread([folder file],'time'); %days since 1950-01-01 00:00:00
time=time+datenum(1950,01,01);
nt=length(time);


% Dimension is different for RACMO
Uwind=ncread([folder file],'u10m',[idr(1) jdr(1) 1 1],[irl jrl  1 Inf]);

% Set new time-limit
time_lim=[datenum(2016,01,01) datenum(2017,01,01)];
ind1=find(time>=time_lim(1));
ind2=find(time<=time_lim(2));
ind=intersect(ind1,ind2);

time=time(ind);check_time=datestr(time);
Uwind=Uwind(:,:,1,ind); %Dimension is different for RACMO
nt=length(time);

UW=NaN(nc,nt);

for i=1:nt
      
    TT=squeeze(Uwind(:,:,i)); %202 x 216 x (each time step=1)
    TT=TT(:); % 202 * 216 = 43632 x 1 
    UW(:,i)=sum(TT(cindex).*ccoef,2);
    clear TT
end 

figure(1);clf
lon=ncread('u10m.nc','lon');lat=ncread('u10m.nc','lat');
lon=lon(idr,jdr);lat=lat(idr,jdr);
%Plot cropped long-lat and data
pcolor(lon,lat,Uwind(:,:,1,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M;[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
scatter(lonc,latc,10,UW(:,1))
colorbar;colormap jet;
ylim([73 87]);


save(['./Racmo_Forc/2016/UW_RAC_' num2str(year)],'time','UW','-v7.3')
clear Uwind UW


%% V-component (in m/s)
file='v10m.nc';
Vwind=ncread([folder file],'v10m',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);
Vwind=Vwind(:,:,ind);

VW=NaN(nc,nt);
 
for i=1:nt
    TT=squeeze(Vwind(:,:,i));
    TT=TT(:);
    VW(:,i)=sum(TT(cindex).*ccoef,2);
    clear TT    
end

figure(2);clf
%Plot cropped long-lat and data
pcolor(lon,lat,Vwind(:,:,1,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M;[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
scatter(lonc,latc,10,VW(:,1))
colorbar;colormap jet;
ylim([73 87]);


save(['./Racmo_Forc/2016/VW_RAC_' num2str(year)],'time','VW','-v7.3')
clear Vwind VW
