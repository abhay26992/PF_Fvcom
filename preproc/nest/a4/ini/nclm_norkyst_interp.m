
function[]=nclm_norkyst_interp(romspath,ncname,year,month)


%%%%%%%%path and files%%%%%%%%%%%
%romspath='/global/work/apn/A4_modelruns/A4_nudging/linked_with_dates/';
%ncname='a4_avg_';
%year=['2014';'2015';'2016']; %'start_year';...;'end_year'
%month=['06';'07';'08'];

%USAGE:
%nclm_norkyst_interp(romspath,ncname,year,month);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%nclm_norkyst_interp.m interpolates roms results onto fvcom node and cells, daily

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





%%%%%start from 'the specified year and months'%%%%%%%%%%%%%%

for i=1
  for j=1:size(month,1) %To allow the code to scan through the '3' AOTIM months
         
        % The rest of the data
        A1=dir([romspath  ncname  year(i,:) '-' month(j,:) '-' '*']);
        
        nf=length(A1);
        
        %nclm.zeta=NaN(0,0);
        %nclm.u=NaN(0,0,0);
        %nclm.v=NaN(0,0,0);
        %nclm.ubar=NaN(0,0);
        %nclm.vbar=NaN(0,0);
        %nclm.salt=NaN(0,0,0);
        %nclm.temp=NaN(0,0,0);
        
        for k=1:nf

            %%%%%%%
            romsfile=[romspath A1(k).name];
            disp([romsfile])
            
            %zeta
            
            zeta=ncread(romsfile,'zeta',[idr(1) jdr(1) 1],[irl jrl Inf]);
            [dim1,dim2,dim_time]=size(zeta);

            if k==1
                start_index=dim_time;
                time_length=1;
                zeta=zeta(2:end-1,2:end-1,end);
            elseif k==nf
                start_index=1;
                time_length=1;
                zeta=zeta(2:end-1,2:end-1,1);
            else
                start_index=1;
                time_length=1;
                zeta=zeta(2:end-1,2:end-1,1);
            end
            
            %salt
            salt=ncread(romsfile,'salt');
            salt=salt(2:end-1,2:end-1,:,:);
            
            %temp
            temp=ncread(romsfile,'temp');
            temp=temp(2:end-1,2:end-1,:,:);
            
            %u
            u=ncread(romsfile,'u');
            u=0.5*(u(1:end-1,:,:,:)+u(2:end,:,:,:));
            u=u(:,2:end-1,:,:);
            
            %ubar
            ubar=ncread(romsfile,'ubar');
            ubar=0.5*(ubar(1:end-1,:,:)+ubar(2:end,:,:));
            ubar=ubar(:,2:end-1,:);
            
            
            %v
            v=ncread(romsfile,'v');
            v=0.5*(v(:,1:end-1,:,:)+v(:,2:end,:,:));
            v=v(2:end-1,:,:,:);
            
            
            %vbar
            vbar=ncread(romsfile,'vbar');
            vbar=0.5*(vbar(:,1:end-1,:)+vbar(:,2:end,:));
            vbar=vbar(2:end-1,:,:);
            
            
            
            [nx,ny,nz,nt]=size(u);
            
            %%%%%%%%
            
            ZETA=NaN(nn,nt);
            UBAR=NaN(nc,nt);
            VBAR=NaN(nc,nt);
            SALT=NaN(nn,nz,nt);
            TEMP=NaN(nn,nz,nt);
            U=NaN(nc,nz,nt);
            V=NaN(nc,nz,nt);
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %disp(['horizontal interpolation'] )
            disp([year(i,:) ' ' month(j,:) ' ' num2str(k) ' of ' num2str(nf)])             
           
            for n=1:nt  %each time interval
                
                %zeta 
                Z0=squeeze(zeta(:,:,n));
                Z0=Z0(ocn_r); % only choose the points in water
                ZETA(:,n)=sum(ncoef_r.*Z0(nind_r),2);
                      
                %ubar
                UB=squeeze(ubar(:,:,n));
                
                UB=UB(ocn_r); % only choose the ocean points
                UBAR(:,n)=sum(ncoef_r2c.*UB(nind_r2c),2);
                
                
                %vbar
                
                VB=squeeze(vbar(:,:,n));
                
                VB=VB(ocn_r); % only choose the ocean points
                VBAR(:,n)=sum(ncoef_r2c.*VB(nind_r2c),2);
                
                %%% convert to ubar and vbar to true north and east direction
                
                [UBARGEO(:,n),VBARGEO(:,n)]=vel_xy2geo(UBAR(:,n),VBAR(:,n),ANGLE);
                
                %u,v
                UU=NaN(nc,nz);
                VV=NaN(nc,nz);
                UUGEO=NaN(nc,nz);
                VVGEO=NaN(nc,nz);
                SS=NaN(nn,nz);
                TT=NaN(nn,nz);
                
                
                for l=1:nz
                    
                    S0=squeeze(salt(:,:,l,n));
                    S0=S0(ocn_r);
                    SS(:,l)=sum(ncoef_r.*S0(nind_r),2);
                    
                    
                    T0=squeeze(temp(:,:,l,n));
                    T0=T0(ocn_r);
                    TT(:,l)=sum(ncoef_r.*T0(nind_r),2);
                    
                    
                    U0=squeeze(u(:,:,l,n));
                    U0=U0(ocn_r);
                    UU(:,l)=sum(ncoef_r2c.*U0(nind_r2c),2);
                    
                    
                    V0=squeeze(v(:,:,l,n));
                    V0=V0(ocn_r);
                    VV(:,l)=sum(ncoef_r2c.*V0(nind_r2c),2);
                    
                    
                    [UUGEO(:,l),VVGEO(:,l)]=vel_xy2geo(UU(:,l),VV(:,l),ANGLE);
                    
                end
                
                SALT(:,:,n)=SS;
                TEMP(:,:,n)=TT;
                U(:,:,n)=UUGEO;
                V(:,:,n)=VVGEO;
            end
            
            
            
            
            
            
            
            clear zeta salt temp SS TT S0 T0  U0 V0 ubar vbar u v UU VV  nx ny
            
            
            
            
            %%%vertical interpolation
            %disp(['vertical interpolation'] )
            
            salt=NaN(nn,kb-1,nt);
            temp=NaN(nn,kb-1,nt);
            u=NaN(nc,kb-1,nt);
            v=NaN(nc,kb-1,nt);
            
            %%%%%%%%%%%%%%%%%rho points
            
            
            for n=1:nn
                
                coef=squeeze(zcoef_n(n,:,:));
                ind=squeeze(zind_n(n,:,:));
                for l=1:nt
                    
                    SS=squeeze(SALT(n,:,l));
                    TT=squeeze(TEMP(n,:,l));
                    
                    salt(n,:,l)=sum(coef.*SS(ind),2);
                    temp(n,:,l)=sum(coef.*TT(ind),2);
                end
            end
            
            
            clear SALT TEMP SS TT coef ind
            
            
            
            %%%%%%%%uv points%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            for n=1:nc
                coef=squeeze(zcoef_c(n,:,:));
                ind=squeeze(zind_c(n,:,:));
                
                
                
                for  l=1:nt
                    
                    UU=squeeze(U(n,:,l));
                    VV=squeeze(V(n,:,l));
                    
                    u(n,:,l)=sum(coef.*UU(ind),2);
                    v(n,:,l)=sum(coef.*VV(ind),2);
                    
                end
                
            end
            
            
            
            clear U V UU VV
            
            if exist('nclm')
                nclm.u=cat(3,nclm.u,u); clear u
                nclm.v=cat(3,nclm.v,v); clear v
                nclm.ubar=cat(2,nclm.ubar,UBARGEO); clear UBARGEO
                nclm.vbar=cat(2,nclm.vbar,VBARGEO); clear VBARGEO
                nclm.salt=cat(3,nclm.salt,salt); clear salt
                nclm.temp=cat(3,nclm.temp,temp); clear temp
                nclm.zeta=cat(2,nclm.zeta,ZETA); clear ZETA
             else
                nclm.u(1:size(u,1),1:size(u,2),1)=u; clear u
                nclm.v(1:size(v,1),1:size(v,2),1)=v; clear v
                nclm.salt(1:size(salt,1),1:size(salt,2),1)=salt; clear salt
                nclm.temp(1:size(temp,1),1:size(temp,2),1)=temp; clear temp
                nclm.ubar(1:size(UBARGEO,1),1)=UBARGEO; clear UBARGEO
                nclm.vbar(1:size(VBARGEO,1),1)=VBARGEO; clear VBARGEO
                nclm.zeta(1:size(ZETA,1),1)=ZETA; clear ZETA
             end
           
        end
        
        %nclm.time=datenum(str2num(year(i,:)),str2num(month(j,:)),1,0:(size(nclm.zeta,2)-1),0,0);
        
        nclm.time=zeros(1,nf);
        for kk=1:nf
            nclm.time(1,kk)=datenum(str2num(year(i,:)),str2num(month(j,:)),kk); %In nf, we store the number of days
        end

        %casename=get_cstr;
        casename='pf_o_1' 
        
        filename=[casename '_nest_',year(i,:),'_',month(j,:),'_all_days'];
        disp([filename])
        
        save(['/home/abhay/Matlab_Repository/Petermann_Bathy/Setup_Nesting/fv_tools/FV_R_PFjord/' filename],'nclm','-v7.3')
        
        clear nclm
        
        
        
    end
end



