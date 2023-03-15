function write_hf_nest(ngrd,Mobj,nest_struct,fileprefix)

use nest_struct

use ngrd  
%Note ngrd -> roms (xc(614x1),yc,xn(618x1),yn,nv(614x3),lonn,latn,lonc,latc,h,fvangle,R,nid,cid,oend1,oend2)
nn=length(xn);
nc=length(xc);
kb=get_kb;
nt=length(hs_time); %hs_time = upsampled hourly time

%
siglay=Mobj.siglay(nid,:); % fv4 
siglev=Mobj.siglev(nid,:); % fv4

%
W=zeros(nn,kb,size(hs_zeta,2));

disp('start dumping')

ncid=netcdf.create([fileprefix, '.nc'], 'clobber');


%% 1. Define dimensions:

timedimid=netcdf.defDim(ncid,'time',netcdf.getConstant('UNLIMITED'));
nodedimid=netcdf.defDim(ncid,'node',nn);
celldimid=netcdf.defDim(ncid,'nele',nc);
threedimid=netcdf.defDim(ncid,'three',3);
levdimid=netcdf.defDim(ncid,'siglev',kb);
laydimid=netcdf.defDim(ncid,'siglay',kb-1);


%% 2. Define variables:

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
uaid=netcdf.defVar(ncid,'ua','float',[celldimid,timedimid]);
vaid=netcdf.defVar(ncid,'va','float',[celldimid,timedimid]);
uid=netcdf.defVar(ncid,'u','float',[celldimid,laydimid,timedimid]);
vid=netcdf.defVar(ncid,'v','float',[celldimid,laydimid,timedimid]);
tempid=netcdf.defVar(ncid,'temp','float',[nodedimid,laydimid,timedimid]);
saltid=netcdf.defVar(ncid,'salinity','float',[nodedimid,laydimid,timedimid]);
hywid=netcdf.defVar(ncid,'hyw','float',[nodedimid,levdimid,timedimid]);

wcid=netcdf.defVar(ncid,'weight_cell','float',[celldimid,timedimid]);
wnid=netcdf.defVar(ncid,'weight_node','float',[nodedimid,timedimid]);
  
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

%% 3. Load data into netcdf file variables
netcdf.putVar(ncid,nvid,nv);
netcdf.putVar(ncid,xnid,x);
netcdf.putVar(ncid,ynid,y);
netcdf.putVar(ncid,xcid,xc);
netcdf.putVar(ncid,ycid,yc);
netcdf.putVar(ncid,lonid,lon);
netcdf.putVar(ncid,latid,lat);
netcdf.putVar(ncid,loncid,lonc);
netcdf.putVar(ncid,latcid,latc);
time=hs_time; %Already in FVCOM time
Itime=floor(time);
Itime2=(time-Itime)*86400*1000;
netcdf.putVar(ncid,timeid,0,length(time),time);
netcdf.putVar(ncid,itimeid,Itime);
netcdf.putVar(ncid,itime2id,Itime2);
netcdf.putVar(ncid,saltid,hs_salinity);
netcdf.putVar(ncid,tempid,hs_temp);
netcdf.putVar(ncid,hywid,hs_hyw);
netcdf.putVar(ncid,zetaid,hs_zeta);
netcdf.putVar(ncid,hid,h);
hc=(h(nv(:,1))+h(nv(:,2))+h(nv(:,3)))/3;%fv4
netcdf.putVar(ncid,hcid,hc);
netcdf.putVar(ncid,siglayid,siglay);
netcdf.putVar(ncid,siglevid,siglev);
slayc=(siglay(nv(:,1),:)+siglay(nv(:,2),:)+siglay(nv(:,3),:))/3; %fv4
slevc=(siglev(nv(:,1),:)+siglev(nv(:,2),:)+siglev(nv(:,3),:))/3;  %fv4
netcdf.putVar(ncid,siglaycid,slayc);
netcdf.putVar(ncid,siglevcid,slevc);
netcdf.putVar(ncid,uid,hs_u);
netcdf.putVar(ncid,vid,hs_v);
netcdf.putVar(ncid,uaid,hs_ua);
netcdf.putVar(ncid,vaid,hs_va);
netcdf.putVar(ncid,wcid,hs_weight_cell);
netcdf.putVar(ncid,wnid,hs_weight_node);
  
netcdf.close(ncid);


