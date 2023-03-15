function[]=cice_forc(cicepath,ncname,year,month)

% Following modifications are introduced:

% Converting it to a function, and one that can day-wise (A4 time step)
% loop over years of data.

% Fixing horrible bugs

%%%%%%%%path and files%%%%%%%%%%%
%cicepath='/global/work/apn/A4_modelruns/A4_nudging/cice/';
%cicepath='C:\PhD\FVCOM\Matlab_Repository\Petermann_Bathy\Setup_Nesting\Roms_a4_daily_nest\FV_R_PFjord\';
%ncname='iceh.';
%year=['2014';'2015';'2016']; %'start_year';...;'end_year'
%month=['06';'07';'08']; %'start_month';...;'end_month'

%USAGE:
%cice_forc(cicepath,ncname,year,month);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%iclm_norkyst_interp.m interpolates roms results onto fvcom node and cells, daily

%clear all
%close all

% Using a var name identical to an inbuilt Matlab fn is horrible practice. 
% Fixed it by assigning/initialising "angle" as empty var.

angle=[]; 
load ROGrd

use ROGrd  % cid angle h latc latn lonc lomn nid nv R cc xn yc yn


% Load index for the subdomain
load rind

use rind
irl=length(idr);
jrl=length(jdr);

iul=length(idu);
jul=length(jdu);

ivl=length(idv);
jvl=length(jdv);

% Load index and coefficient of four nearest points to fvcom node

load nearest4_r  % coef_r ind_r

load nearest4_r2c % coef_r2c ind_r2c

% Mobj
load ./M.mat
xn=Mobj.x;yn=Mobj.y;
xc=Mobj.xc;yc=Mobj.yc;
nn=length(xn);
nc=length(xc);

%% Loop :


%%%%%%%%%%%%%%%  start from 'the specified year and months'  %%%%%%%%%%%%%%

for i=1
  for j=1:size(month,1) %To allow the code to scan through the months
         
        % The rest of the data
        A1=dir([cicepath  ncname  year(i,:) '-' month(j,:) '-' '*']);
        
        nf=length(A1); %Number of days in that month
        
        %nclm.zeta=NaN(0,0);
        %nclm.u=NaN(0,0,0);
        %nclm.v=NaN(0,0,0);
        %nclm.ubar=NaN(0,0);
        %nclm.vbar=NaN(0,0);
        %nclm.salt=NaN(0,0,0);
        %nclm.temp=NaN(0,0,0);
        
        for k=1:nf   % Scan through all days of month, one day at a time

           
            romsfile=[cicepath A1(k).name];
            disp([romsfile])
            
            
            %aicen
            aicen_d=ncread(romsfile,'aicen_d',[idr(1) jdr(1) 1 1],[irl jrl Inf Inf]);
            aicen_d=aicen_d(2:end-1,2:end-1,:,:);

            
            %vicen
            vicen_d=ncread(romsfile,'vicen_d',[idr(1) jdr(1) 1 1],[irl jrl Inf Inf]); %Change -> Using "rind" to crop
            vicen_d=vicen_d(2:end-1,2:end-1,:,:);
            
            %vsnon
            vsnon_d=ncread(romsfile,'vsnon_d',[idr(1) jdr(1) 1 1],[irl jrl Inf Inf]);
            vsnon_d=vsnon_d(2:end-1,2:end-1,:,:); % 2:end
            vsnon_r=vsnon_d;
            
            %sice_d
            sice_d=ncread(romsfile,'sice_d',[idr(1) jdr(1) 1],[irl jrl Inf]);
            sice_d=sice_d(2:end-1,2:end-1,:); % 2:end
            
            %tsfc_d
            Tsfc_d=ncread(romsfile,'Tsfc_d',[idr(1) jdr(1) 1],[irl jrl Inf]);
            Tsfc_d=Tsfc_d(2:end-1,2:end-1,:); % 2:end

            %uvel
            uvel_d=ncread(romsfile,'uvel_d',[idu(1) jdu(1) 1],[iul jul Inf]);
            uvel_d=0.5*(uvel_d(1:end-1,:,:)+uvel_d(2:end,:,:));
            uvel_d=uvel_d(:,2:end-1,:);

            %vvel
            vvel_d=ncread(romsfile,'vvel_d',[idv(1) jdv(1) 1],[ivl jvl Inf]);
            vvel_d=0.5*(vvel_d(:,1:end-1,:)+vvel_d(:,2:end,:));
            vvel_d=vvel_d(2:end-1,:,:);
            
            %%%%%%%%%%%%%%% vel_xy2geo is done here %%%%%%%%%%
            [uvelgeo,vvelgeo]=vel_xy2geo(uvel_d,vvel_d,angle);

            %tlon
            TLON=ncread(romsfile,'TLON',[idr(1) jdr(1)],[irl jrl]);
            TLON=TLON(2:end-1,2:end-1); % 2:end
            
            %tlat
            TLAT=ncread(romsfile,'TLAT',[idr(1) jdr(1)],[irl jrl]);
            TLAT=TLAT(2:end-1,2:end-1); % 2:end
            
            % Plot:
            
            %%% Plots:
            %keyboard
            
            %load ./M.mat;
            %figure(1);clf
            %plot(TLON-360,TLAT);hold on;scatter(Mobj.lon,Mobj.lat);
            
            
            %figure(2);clf;pcolor(uvel_d(:,:));colorbar;
            
            %figure(3);clf;load ./M.mat;
            %pcolor(TLON-360,TLAT,Tsfc_d);colorbar;colormap jet;hold on;
            %title('snow/ice surface temp in Celsius');
            %scatter(Mobj.lon,Mobj.lat,'o','k');

            
            [nx,ny,nz,nt]=size(aicen_d);
            
            %%%%%%%%
            
            % Dimension (lat,lon,time) -> Tsfc_d sice_d uvel_d and vvel_d
            TSFC_D=NaN(nn,nt);
            SICE_D=NaN(nn,nt);
            UVEL_D=NaN(nc,nt);
            VVEL_D=NaN(nc,nt);
            
            % Dimension (lat,lon,thickness,time) -> aicen_d vicen_d vsnon_d
            AICEN_D=NaN(nn,nz,nt);
            VICEN_D=NaN(nn,nz,nt);
            VSNON_D=NaN(nn,nz,nt);
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp(['horizontal interpolation'] )
            disp([year(i,:) ' ' month(j,:) ' ' num2str(k) ' of ' num2str(nf)])             
           
            for n=1:nt  %each time interval
                
                %zeta 
                T0=squeeze(Tsfc_d(:,:,n));
                T0=T0(ocn_r); % only choose the points in water
                TSFC_D(:,n)=sum(coef_r.*T0(ind_r),2);
                
                %keyboard
                %scatter(lonn,latn,10,TSFC_D);colorbar;colormap jet;
                %title('snow/ice surface temp in Celsius : Hinterp');
                %load('M.mat');hold on;scatter(Mobj.lon,Mobj.lat,1,'k');
                %tsfc_nan=find(isnan(TSFC_D));
                
                %sice
                SI=squeeze(sice_d(:,:,n));
                SI=SI(ocn_r); % only choose the points in water
                SICE_D(:,n)=sum(coef_r.*SI(ind_r),2);
                %sice_nan=find(isnan(SICE_D));
                      
                %ubar
                UV=squeeze(uvelgeo(:,:,n));
                UV=UV(ocn_r); % only choose the ocean points
                UVEL_D(:,n)=sum(coef_r2c.*UV(ind_r2c),2);
                %uvel_nan=find(isnan(UVEL_D));
                
                %vbar
                VV=squeeze(vvelgeo(:,:,n));
                VV=VV(ocn_r); % only choose the ocean points
                VVEL_D(:,n)=sum(coef_r2c.*VV(ind_r2c),2);
                %vvel_nan=find(isnan(VVEL_D));
                
                %%% convert to ubar and vbar to true north and east direction
                
                %[UVELGEO(:,n),VVELGEO(:,n)]=vel_xy2geo(UVEL_D(:,n),VVEL_D(:,n),angle); 
                %uvgeo_nan=find(isnan(UVELGEO));vvgeo_nan=find(isnan(VVELGEO));
                
                %u,v
                AICEN=NaN(nn,nz);
                VICEN=NaN(nn,nz);
                VSNON=NaN(nn,nz);
                
                
                for l=1:nz %for all the depths from 35 (surface) to 1 (bottom) [for CICE - these are the 5 categories]
                    
                    AI=squeeze(aicen_d(:,:,l,n));
                    AI=AI(ocn_r);
                    AICEN(:,l)=sum(coef_r.*AI(ind_r),2);
                    %AICEN_NaN=find(isnan(AICEN));
                    
                    VI=squeeze(vicen_d(:,:,l,n));
                    VI=VI(ocn_r);
                    VICEN(:,l)=sum(coef_r.*VI(ind_r),2);
                    %VICEN_NaN=find(isnan(VICEN));
                    
                    VS=squeeze(vsnon_d(:,:,l,n));
                    VS=VS(ocn_r);
                    VSNON(:,l)=sum(coef_r.*VS(ind_r),2);
                    %VSNON_NaN=find(isnan(VSNON));
                   
                end
                
                AICEN_D(:,:,n)=AICEN;
                VICEN_D(:,:,n)=VICEN;
                VSNON_D(:,:,n)=VSNON;
            end
            
            %End of horizontal interp -> for each time interval and for all
            %the depth. So, for the given day (nc file), and all the
            %depths.
            
            clear aicen_d vicen_d vsnon_d Tsfc_d uvel_d vvel_d sice_d AI VI VS AICEN VICEN VSNON nx ny
            
            %pick a category, or ...
            %aicen_d=AICEN_D(:,5);vicen_d=VICEN_D(:,5);vsnon_d=VSNON_D(:,5);
            
            % use all categories ...
            aicen_d=AICEN_D(:,:);vicen_d=VICEN_D(:,:);vsnon_d=VSNON_D(:,:);
            
            
%% Skip vertical interp for CICE   

            %disp(['vertical interpolation'] )
            
%             aicen_d=NaN(nn,kb-1,nt);
%             vicen_d=NaN(nn,kb-1,nt);
%             vsnon_d=NaN(nn,kb-1,nt);
%             
%             %rho points
%             
%             
%             for n=1:nn %One node point (and all depth levels) at a time
%                 
%                 coef=squeeze(zcoef_n(n,:,:)); %zcoef_n -> 618 x 23 x 2      coef -> (23 x 2) for n = 1:618 
%                 %[For each day, we go from 1 : 618 using 23 x 2 coef
%                 %values. Each of the 618 points have 23 x 2 coef values.
%                 %All looks good for day 1, 2 ... warming for day 20]
%                 ind=squeeze(zind_n(n,:,:));
%                 for l=1:nt
%                     
%                     AICEN=squeeze(AICEN_D(n,:,l));  
%                     VICEN=squeeze(VICEN_D(n,:,l)); % node n = 1:618, all depth, l=1 (always)
%                     VSNON=squeeze(VSNON_D(n,:,l));
%                     
%                     aicen_d(n,:,l)=sum(coef.*AICEN(ind),2);
%                     vicen_d(n,:,l)=sum(coef.*VICEN(ind),2);
%                     vsnon_d(n,:,l)=sum(coef.*VSNON(ind),2); %Filling up temp one node point (and all depth levels) at a time
%                 end
%             end
%                  
%             %Insert Breakpoint here - Debug
%             clear AICEN VICEN VSNON AICEN_D VICEN_D VSNON_D coef ind

%% Populate the struct:
            
            if exist('iclm')
                iclm.u=cat(2,iclm.u,UVEL_D); clear UVELGEO
                iclm.v=cat(2,iclm.v,VVEL_D); clear VVELGEO
                iclm.tsfc=cat(2,iclm.tsfc,TSFC_D); clear TSFC_D
                iclm.sice=cat(2,iclm.sice,SICE_D); clear SICE_D
                iclm.aicen=cat(3,iclm.aicen,aicen_d); clear aicen_d
                iclm.vicen=cat(3,iclm.vicen,vicen_d); clear vicen_d
                iclm.vsnon=cat(3,iclm.vsnon,vsnon_d); clear vsnon_d
            else
                iclm.aicen(1:size(aicen_d,1),1:size(aicen_d,2),1)=aicen_d; clear aicen_d
                iclm.vicen(1:size(vicen_d,1),1:size(vicen_d,2),1)=vicen_d; clear vicen_d
                iclm.vsnon(1:size(vsnon_d,1),1:size(vsnon_d,2),1)=vsnon_d; clear vsnon_d
                iclm.u(1:size(UVEL_D,1),1)=UVEL_D; clear UVELGEO
                iclm.v(1:size(VVEL_D,1),1)=VVEL_D; clear VVELGEO
                iclm.tsfc(1:size(TSFC_D,1),1)=TSFC_D; clear TSFC_D
                iclm.sice(1:size(SICE_D,1),1)=SICE_D; clear SICE_D
             end
           
        end    %"nf" for loop ends here -> having generated data for one full month (nf number of days)
        
        %iclm.u(:)=0;iclm.v(:)=0;
        %iclm.ubar(:)=0;iclm.vbar(:)=0;
        %nclm.time=datenum(str2num(year(i,:)),str2num(month(j,:)),1,0:(size(nclm.zeta,2)-1),0,0);
        
        iclm.time=zeros(1,nf);
        for kk=1:nf
            iclm.time(1,kk)=datenum(str2num(year(i,:)),str2num(month(j,:)),kk); %In nf, we store the number of days
        end

        %casename=get_cstr;
        casename='pf_o_1' 
        
        filename=[casename '_icenudge_',year(i,:),'_',month(j,:),'_all_days'];
        disp([filename])
        
        save(['C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/' filename],'iclm','-v7.3')
        
        clear iclm
        
    end
end