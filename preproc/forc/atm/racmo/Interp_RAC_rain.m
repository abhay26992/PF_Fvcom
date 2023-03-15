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
  
%% Read in racmo data, pay attention to the signs and units%%%%%%%%%%%

%% Precipitation:

folder='C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Racmo_atm_forc/RACMO_Forcing/';
file='precip.nc';

Rain=ncread([folder file],'precip',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);

% 1 kg of precipitation over 1 m2 area has 1 mm thickness ( 1 kg m-2 s-1 =
% 1 mm/s) -> change unit from kg m-2 s-1 (mm/s) to m/s 

%
%Rain=Rain.*3*60*60; % We have 3 hourly averaged values -> Going from kg m-2 s-1 (mm/s) to mm w.e. (or total rainfall in 3 hours)
%Rain=Rain.*1e-3; % Go from mm to m 
%Rain=Rain./(3*60*60); % Convert to m/s 

Rain=Rain.*1e-3; % Go from mm/s to m/s 


% Change time for RACMO (Days since "1950-01-01")
time=ncread([folder file],'time'); %days since 1950-01-01 00:00:00
time=time+datenum(1950,01,01);

% Set new time-limit
time_lim=[datenum(2016,01,01) datenum(2017,01,01)];
ind1=find(time>=time_lim(1));
ind2=find(time<=time_lim(2));
ind=intersect(ind1,ind2);

time=time(ind);check_time=datestr(time);
Rain=Rain(:,:,1,ind); %Dimension is different for RACMO
nt=length(time);


%%%%%% Set the "constant" negative [-7.2338e-10] rainfall values to 0 %%%%%
Rain(Rain<0)=0; 
check=Rain(Rain<0); % check=Rain(isnan(Rain)); check=~nnz(Rain);


RA=NaN(nn,nt);
 for i=1:nt
         
     TT=squeeze(Rain(:,:,i));
     TT=TT(:);
     RA(:,i)=sum(TT(nindex).*ncoef,2);
     clear TT
         
 end
 
figure(1);clf
lon=ncread('precip.nc','lon');lat=ncread('precip.nc','lat');
lon=lon(idr,jdr);lat=lat(idr,jdr);
%Plot cropped long-lat and data
pcolor(lon,lat,Rain(:,:,1,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,10,RA(:,1))
colorbar;colormap jet;
ylim([73 87]);
 
save(['./Racmo_Forc/2016/RA_RAC_' num2str(year)],'RA','time','-v7.3') 
 
 
 %% Evaporation:
file='evap.nc'; 
Ev=ncread([folder file],'evap',[idr(1) jdr(1) 1 1],[irl jrl 1 Inf]);
Ev=Ev(:,:,1,ind);
Ev=Ev.*1e-3;

for i=1:nt
    
    TT=squeeze(Ev(:,:,i));
    TT=TT(:);
    Evap(:,i)=sum(TT(nindex).*ncoef,2);
    clear TT
    
end

figure(2);clf
%Plot cropped long-lat and data
pcolor(lon,lat,Ev(:,:,1,1));shading flat;hold on; 
%Overlay Interpolated result - FVCOM TA
load ./M
scatter(Mobj.lon,Mobj.lat,10,Evap(:,1))
colorbar;colormap jet;
ylim([73 87]); 

save(['./Racmo_Forc/2016/Evap_RAC_' num2str(year)],'time','Evap','-v7.3')