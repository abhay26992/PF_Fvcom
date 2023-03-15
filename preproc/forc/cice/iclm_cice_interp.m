function[]=iclm_cice_interp(cicepath,ncname,year,month)


%%%%%%%%path and files%%%%%%%%%%%
%cicepath='/global/work/apn/A4_modelruns/A4_nudging/cice/';
%cicepath='C:\PhD\FVCOM\Matlab_Repository\Petermann_Bathy\Setup_Nesting\Roms_a4_daily_nest\FV_R_PFjord\';
%ncname='iceh.';
%year=['2014';'2015';'2016']; %'start_year';...;'end_year'
%month=['06';'07';'08'];

%USAGE:
%iclm_cice_interp(cicepath,ncname,year,month);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%iclm_norkyst_interp.m interpolates roms results onto fvcom node and cells, daily

%clear all
%close all


load ROGrd

use ROGrd  % cid angle h latc latn lonc lomn nid nv R cc xn yc yn


%%%%%%%load index for the subdomain
load rind

use rind
irl=length(idr);
jrl=length(jdr);

iul=length(idu);
jul=length(jdu);

ivl=length(idv);
jvl=length(jdv);



%%%%%%%%%%

load ngrd
use ngrd

%%%%%load index and coefficient of four nearest points to fvcom node

load nearest4_r  % coef_r ind_r

load Vinterp    % ZN ZC zind_c zind_n
use Vinterp

load nearest4_r2c

%%%%load nudging zone point


kb=get_kb;



nn=length(nid);
nc=length(cid);



nind_r=ind_r(nid,:);

ncoef_r=coef_r(nid,:);


nind_r2c=ind_r2c(cid,:);
ncoef_r2c=coef_r2c(cid,:);


zind_n=zind_n(nid,:,:);
zcoef_n=zcoef_n(nid,:,:);

zind_c=zind_c(cid,:,:);
zcoef_c=zcoef_c(cid,:,:);

%%%interpolate angle

ANGLE=sum(ncoef_r2c.*angler(nind_r2c),2);





%%%%%%%%%%%%%%%  start from 'the specified year and months'  %%%%%%%%%%%%%%

for i=1
  for j=1:size(month,1) %To allow the code to scan through the '3' AOTIM months
         
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

            %%%%%%%
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
            uvel_d=ncread(romsfile,'uvel_d',[idu(1) jdu(1) 1],[iul jul 1]);
            uvel_d=0.5*(uvel_d(1:end-1,:,:)+uvel_d(2:end,:,:));
            uvel_d=uvel_d(:,2:end-1,:);

            %vvel
            vvel_d=ncread(romsfile,'vvel_d',[idv(1) jdv(1) 1],[ivl jvl 1]);
            vvel_d=0.5*(vvel_d(:,1:end-1,:)+vvel_d(:,2:end,:));
            vvel_d=vvel_d(2:end-1,:,:);

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
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp(['horizontal interpolation'] )
            disp([year(i,:) ' ' month(j,:) ' ' num2str(k) ' of ' num2str(nf)])             
           
            for n=1:nt  %each time interval
                
                %zeta 
                T0=squeeze(Tsfc_d(:,:,n));
                T0=T0(ocn_r); % only choose the points in water
                TSFC_D(:,n)=sum(ncoef_r.*T0(nind_r),2);
                
                %keyboard
                %scatter(lonn,latn,10,TSFC_D);colorbar;colormap jet;
                %title('snow/ice surface temp in Celsius : Hinterp');
                %load('M.mat');hold on;scatter(Mobj.lon,Mobj.lat,1,'k');
                %tsfc_nan=find(isnan(TSFC_D));
                
                %sice
                SI=squeeze(sice_d(:,:,n));
                SI=SI(ocn_r); % only choose the points in water
                SICE_D(:,n)=sum(ncoef_r.*SI(nind_r),2);
                %sice_nan=find(isnan(SICE_D));
                      
                %ubar
                UV=squeeze(uvel_d(:,:,n));
                UV=UV(ocn_r); % only choose the ocean points
                UVEL_D(:,n)=sum(ncoef_r2c.*UV(nind_r2c),2);
                %uvel_nan=find(isnan(UVEL_D));
                
                %vbar
                VV=squeeze(vvel_d(:,:,n));
                VV=VV(ocn_r); % only choose the ocean points
                VVEL_D(:,n)=sum(ncoef_r2c.*VV(nind_r2c),2);
                %vvel_nan=find(isnan(VVEL_D));
                
                %%% convert to ubar and vbar to true north and east direction
                
                [UVELGEO(:,n),VVELGEO(:,n)]=vel_xy2geo(UVEL_D(:,n),VVEL_D(:,n),ANGLE); %Using ANGLE (614 x 1), instead using angle (324x300)?
                %uvgeo_nan=find(isnan(UVELGEO));vvgeo_nan=find(isnan(VVELGEO));
                
                %u,v
                AICEN=NaN(nn,nz);
                VICEN=NaN(nn,nz);
                VSNON=NaN(nn,nz);
                
                
                for l=1:nz %for all the depths from 35 (surface) to 1 (bottom)
                    
                    AI=squeeze(aicen_d(:,:,l,n));
                    AI=AI(ocn_r);
                    AICEN(:,l)=sum(ncoef_r.*AI(nind_r),2);
                    %AICEN_NaN=find(isnan(AICEN));
                    
                    VI=squeeze(vicen_d(:,:,l,n));
                    VI=VI(ocn_r);
                    VICEN(:,l)=sum(ncoef_r.*VI(nind_r),2);
                    %VICEN_NaN=find(isnan(VICEN));
                    
                    VS=squeeze(vsnon_d(:,:,l,n));
                    VS=VS(ocn_r);
                    VSNON(:,l)=sum(ncoef_r.*VS(nind_r),2);
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
            
            %pick a category
            %aicen_d=AICEN_D(:,5);vicen_d=VICEN_D(:,5);vsnon_d=VSNON_D(:,5);
            
            aicen_d=AICEN_D(:,:);vicen_d=VICEN_D(:,:);vsnon_d=VSNON_D(:,:);
            
            
            %%% vertical interpolation
            %disp(['vertical interpolation'] )
            
%             aicen_d=NaN(nn,kb-1,nt);
%             vicen_d=NaN(nn,kb-1,nt);
%             vsnon_d=NaN(nn,kb-1,nt);
%             
%             %%%%%%%%%%%%%%%%%rho points
%             
%             
%             for n=1:nn %One node point (and all depth levels) at a time
%                 
%                 coef=squeeze(zcoef_n(n,:,:)); %zcoef_n -> 618 x 23 x 2      coef -> (23 x 2) for n = 1:618 
%                 %[For each day, we go from 1 : 618 using 23 x 2 coef
%                 %values. Each of the 618 points have 23 x 2 coef values.
%                 %All good for day 1, 2 ... warming for day 20]
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
%             
            if exist('iclm')
                iclm.u=cat(2,iclm.u,UVELGEO); clear UVELGEO
                iclm.v=cat(2,iclm.v,VVELGEO); clear VVELGEO
                iclm.tsfc=cat(2,iclm.tsfc,TSFC_D); clear TSFC_D
                iclm.sice=cat(2,iclm.sice,SICE_D); clear SICE_D
                iclm.aicen=cat(3,iclm.aicen,aicen_d); clear aicen_d
                iclm.vicen=cat(3,iclm.vicen,vicen_d); clear vicen_d
                iclm.vsnon=cat(3,iclm.vsnon,vsnon_d); clear vsnon_d
            else
                iclm.aicen(1:size(aicen_d,1),1:size(aicen_d,2),1)=aicen_d; clear aicen_d
                iclm.vicen(1:size(vicen_d,1),1:size(vicen_d,2),1)=vicen_d; clear vicen_d
                iclm.vsnon(1:size(vsnon_d,1),1:size(vsnon_d,2),1)=vsnon_d; clear vsnon_d
                iclm.u(1:size(UVELGEO,1),1)=UVELGEO; clear UVELGEO
                iclm.v(1:size(VVELGEO,1),1)=VVELGEO; clear VVELGEO
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
        
        filename=[casename '_icenest_',year(i,:),'_',month(j,:),'_all_days'];
        disp([filename])
        
        save(['C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Setup_Nesting/Roms_a4_daily_nest' filename],'iclm','-v7.3')
        
        clear iclm
        
        
        
    end
end