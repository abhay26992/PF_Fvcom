function write_FVCOM_nest(ngrd, fileprefix, nclm,time,nest_type,fclm)
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
if nargin==6
  fabm=true;
end

use ngrd
nn=length(xn);
nc=length(xc);
kb=get_kb;
nt=length(time);


use nclm
if exist('ua','var')&~exist('ubar','var')
 ubar=ua;
 vbar=va;
end

%w
W=zeros(nn,kb,nt);




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
time=time-datenum(1858,11,17,0,0,0);
Itime=floor(time);
Itime2=(time-Itime)*86400*1000;
netcdf.putVar(ncid,timeid,0,length(time),time);
netcdf.putVar(ncid,itimeid,Itime);
netcdf.putVar(ncid,itime2id,Itime2);
netcdf.putVar(ncid,saltid,salt);
netcdf.putVar(ncid,tempid,temp);
netcdf.putVar(ncid,uid,u);
netcdf.putVar(ncid,vid,v);
netcdf.putVar(ncid,uaid,ubar);
netcdf.putVar(ncid,vaid,vbar);
netcdf.putVar(ncid,hywid,W);
netcdf.putVar(ncid,zetaid,zeta);
%
if nest_type==3
  w1=2.5e-4;
  w2=2.5e-5;

  [weight_cell,weight_node]=CalcWeights(time,w1,w2);

  %%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%
  load nbe

  weight_node(isonb(nid)==2,:)=1;

  load M
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

