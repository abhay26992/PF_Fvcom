function write_FVCOM_nest(ngrd, fileprefix, nclm,time,nest_type,geo2xy,fclm)
% Write nest climatology out to FVCOM Netcdf  nest file.
%
% write_FVCOM_nest(ngrd, fileprefix,nclm)
%
% DESCRIPTION:
%   Takes the given interpolated climatology data and writes out
%   to a NetCDF file.
%
% INPUT:
%   ngrd - grid information for nesting zone
%   fileprefix - Output NetCDF file prefix (plus path) 
%   nclm - Struct of the data to be written out.
%   time- time of the data
%  

%%%%%%%%%%%%%%get dimensions %%%%%%%%%%%%%%%%%%%%


fabm=false;
if nargin==7
  fabm=true;
end

use ngrd  
%Note ngrd -> roms (xc(614x1),yc,xn(618x1),yn,nv(614x3),lonn,latn,lonc,latc,h,fvangle,R,nid,cid,oend1,oend2)
nn=length(xn);
nc=length(xc);
kb=get_kb;
nt=length(time); %time = nclm.time (roms) -> we have nclm.time starting from June 1,2016

%
load ./M.mat
siglay=Mobj.siglay(nid,:); % fv4 
siglev=Mobj.siglev(nid,:); % fv4
%

use nclm
if exist('ua','var')&~exist('ubar','var') %ua,va = fvcom
 ubar=ua;
 vbar=va;
end

%w
%W=zeros(nn,kb,nt);
W=zeros(nn,kb,size(zeta,2));




%%%%%%%%%%%%%%%%%dumping data to netcdf file%%%%%%%%%%%%%%

disp('start dumping')

ncid=netcdf.create([fileprefix, '.nc'], 'clobber');

%Define dimensions
timedimid=netcdf.defDim(ncid,'time',netcdf.getConstant('UNLIMITED'));
nodedimid=netcdf.defDim(ncid,'node',nn);
celldimid=netcdf.defDim(ncid,'nele',nc);
threedimid=netcdf.defDim(ncid,'three',3);
levdimid=netcdf.defDim(ncid,'siglev',kb);
laydimid=netcdf.defDim(ncid,'siglay',kb-1);
%Define variables
nvid=netcdf.defVar(ncid,'nv','int',[celldimid,threedimid]);
xnid=netcdf.defVar(ncid,'x','float',nodedimid);
ynid=netcdf.defVar(ncid,'y','float',nodedimid);
xcid=netcdf.defVar(ncid,'xc','float',celldimid);
ycid=netcdf.defVar(ncid,'yc','float',celldimid);
lonid=netcdf.defVar(ncid,'lon','float',nodedimid);
latid=netcdf.defVar(ncid,'lat','float',nodedimid);
loncid=netcdf.defVar(ncid,'lonc','float',celldimid);
latcid=netcdf.defVar(ncid,'latc','float',celldimid);
zetaid=netcdf.defVar(ncid,'zeta','float',[nodedimid,timedimid]);
hid=netcdf.defVar(ncid,'h','float',nodedimid);
hcid=netcdf.defVar(ncid,'h_center','float',celldimid); %fv4
siglayid=netcdf.defVar(ncid,'siglay','float',[nodedimid,laydimid]);%fv4
siglevid=netcdf.defVar(ncid,'siglev','float',[nodedimid,levdimid]);%fv4
siglaycid=netcdf.defVar(ncid,'siglay_center','float',[celldimid,laydimid]);%fv4
siglevcid=netcdf.defVar(ncid,'siglev_center','float',[celldimid,levdimid]);%fv4
%zetaid=netcdf.defVar(ncid,'zeta','float',[timedimid,nodedimid]);
uaid=netcdf.defVar(ncid,'ua','float',[celldimid,timedimid]);
vaid=netcdf.defVar(ncid,'va','float',[celldimid,timedimid]);
uid=netcdf.defVar(ncid,'u','float',[celldimid,laydimid,timedimid]);
vid=netcdf.defVar(ncid,'v','float',[celldimid,laydimid,timedimid]);
tempid=netcdf.defVar(ncid,'temp','float',[nodedimid,laydimid,timedimid]);
saltid=netcdf.defVar(ncid,'salinity','float',[nodedimid,laydimid,timedimid]);
hywid=netcdf.defVar(ncid,'hyw','float',[nodedimid,levdimid,timedimid]);
if nest_type==3
  wcid=netcdf.defVar(ncid,'weight_cell','float',[celldimid,timedimid]);
  wnid=netcdf.defVar(ncid,'weight_node','float',[nodedimid,timedimid]);
end
if fabm
  fstr='float';
  for ifabm=1:length(fclm)
    eval(['fabmid' num2str(ifabm) '=netcdf.defVar(ncid,fclm(ifabm).varname,fstr,[nodedimid,laydimid,timedimid]);'])
  end
end
%Time
timeunit='days since 1858-11-17 00:00:00';
timeunit2='msec since 00:00:00';
timeformat='modified julian day (MJD)';
time_zone='UTC';
timeid=netcdf.defVar(ncid,'time','float',timedimid);
netcdf.putAtt(ncid,timeid,'long_name','time')
netcdf.putAtt(ncid,timeid,'units',timeunit)
netcdf.putAtt(ncid,timeid,'format',timeformat)
netcdf.putAtt(ncid,timeid,'time_zone',time_zone)
itimeid=netcdf.defVar(ncid,'Itime','int',timedimid);
netcdf.putAtt(ncid,itimeid,'units',timeunit)
netcdf.putAtt(ncid,itimeid,'format',timeformat)
netcdf.putAtt(ncid,itimeid,'time_zone',time_zone)
itime2id=netcdf.defVar(ncid,'Itime2','int',timedimid);
netcdf.putAtt(ncid,itime2id,'units',timeunit2)
netcdf.putAtt(ncid,itime2id,'time_zone',time_zone)
netcdf.endDef(ncid);
%Load data into netcdf file variables
netcdf.putVar(ncid,nvid,nv);
netcdf.putVar(ncid,xnid,xn);
netcdf.putVar(ncid,ynid,yn);
netcdf.putVar(ncid,xcid,xc);
netcdf.putVar(ncid,ycid,yc);
netcdf.putVar(ncid,lonid,lonn);
netcdf.putVar(ncid,latid,latn);
netcdf.putVar(ncid,loncid,lonc);
netcdf.putVar(ncid,latcid,latc);
time=time-datenum(1858,11,17,0,0,0); %Going back to FVCOM time (from 1 June 2016 to 15 July 0157)
Itime=floor(time);
Itime2=(time-Itime)*86400*1000;
netcdf.putVar(ncid,timeid,0,length(time),time);
netcdf.putVar(ncid,itimeid,Itime);
netcdf.putVar(ncid,itime2id,Itime2);
netcdf.putVar(ncid,saltid,salt);
netcdf.putVar(ncid,tempid,temp);
netcdf.putVar(ncid,hywid,W);
netcdf.putVar(ncid,zetaid,zeta);
netcdf.putVar(ncid,hid,h);
hc=(h(nv(:,1))+h(nv(:,2))+h(nv(:,3)))/3;%fv4
netcdf.putVar(ncid,hcid,hc);
netcdf.putVar(ncid,siglayid,siglay);
netcdf.putVar(ncid,siglevid,siglev);
slayc=(siglay(nv(:,1),:)+siglay(nv(:,2),:)+siglay(nv(:,3),:))/3; %fv4
slevc=(siglev(nv(:,1),:)+siglev(nv(:,2),:)+siglev(nv(:,3),:))/3;  %fv4
netcdf.putVar(ncid,siglaycid,slayc);
netcdf.putVar(ncid,siglevcid,slevc);
% compute ubar,vbar if not exist
if ~exist('ubar','var')
    dz1=slevc(:,1:end-1)-slevc(:,2:end);
    for n=1:nt
    for i=1:nc
        ubar(i,n)=sum(u(i,:,n).*dz1(i,:));
        vbar(i,n)=sum(v(i,:,n).*dz1(i,:));
    end
    end
end

if geo2xy==1
    mfvangle=repmat(fvangle,1,kb-1);
    UBAR=NaN(size(ubar));
    VBAR=UBAR;
    U=NaN(size(u));
    V=NaN(size(v));
    for i=1:size(u,3)
        
        [UBAR(:,i),VBAR(:,i)]=vel_geo2xy(ubar(:,i),vbar(:,i),fvangle);
        [U(:,:,i),V(:,:,i)]=vel_geo2xy(u(:,:,i),v(:,:,i),mfvangle);
    end
    
end

%u=U; clear U
%v=V; clear V
%ubar=UBAR; clear UBAR
%vbar=VBAR; clear VBAR
netcdf.putVar(ncid,uid,u);
netcdf.putVar(ncid,vid,v);
netcdf.putVar(ncid,uaid,ubar);
netcdf.putVar(ncid,vaid,vbar);

%
if nest_type==3
  w1=2.5e-4;
  w2=2.5e-5;

  [weight_cell,weight_node]=CalcWeights(time,w1,w2);


  %%%%%%%%%%%%%%%%%%%%%%
  load ./nbe

  weight_node(isonb(nid)==2,:)=1;

  load ./M
  nv=Mobj.tri;

  nc=length(cid);
  for i=1:nc
    if (isonb(nv(cid(i),1))==2|isonb(nv(cid(i),2))==2|isonb(nv(cid(i),3))==2);
    weight_cell(i,:)=1;
    end
  end



  wc=weight_cell(:,1);
  wn=weight_node(:,1);

  qq=find(wc==1);
  nq=length(qq);

  for i=1:nq
   for  j=1:3
    tmp=nv(cid(qq(i)),j);
    zz=find(nid==tmp);
    if wn(zz)~=1
        weight_node(zz,:)=1;
    end
    
   
    
  end
  end




%%%%%%%%%%%%%%%%%%%%%%%%%

  netcdf.putVar(ncid,wcid,weight_cell);
  netcdf.putVar(ncid,wnid,weight_node);

end


if fabm
  for ifabm=1:length(fclm)
    eval(['netcdf.putVar(ncid,fabmid' num2str(ifabm) ',fclm(' num2str(ifabm) ').data);'])
  end
end


netcdf.close(ncid);




end

