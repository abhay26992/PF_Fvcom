function theta=transport_anglec(Mobj)
%This routine calculates the angle between the x-axis and eastward
%direction. The angle is positive if the x-axis is turned counterclockwise
%relative to eastward direction.
%Usage: theta=anglec(nc,xc,yc,xn,yn,proj)
%Input: proj - sets the projection for calculating lat, lon from x ,y. 
%              proj='pst' or 'utm'
%              nc,xc,yc,xn,yn is output from GridSize.m and readgrid.m

nc=Mobj.nElems;
nn=Mobj.nVerts;
xn=Mobj.x;yn=Mobj.y;
xc=Mobj.xc;yc=Mobj.yc;
tri=Mobj.tri;
ids=[ones(size(tri,1),1),tri];
%[nc,nn]=GridSize;
%[xc,yc,xn,yn]=readgrid(nc,nn);
%if strcmp(proj,'utm')
  %utmzone=repmat('33 W',size(yn,1),1);
  %[latn,lonn]=utm2deg(xn,yn,utmzone);
  %utmzone=repmat('33 W',size(yc,1),1);
  %[latc,lonc]=utm2deg(xc,yc,utmzone);
  
  [latn,lonn]=polarstereo_inv(xn,yn);
  [latc,lonc]=polarstereo_inv(xc,yc);
  
  
  
%elseif strcmp(proj,'pst')
%  [latn,lonn]=polarstereo_inv(xn,yn,[],[],78,10);
%  [latc,lonc]=polarstereo_inv(xc,yc,[],[],78,10);
%else
%  disp('proj must equal utm or pst')
%end
%Read gridfile to get cell and node id's%'
casestr=get_cstr;
%fid=fopen(['input/' casestr '_grd.dat']);
%tline=fgetl(fid);
%tline=fgetl(fid);
%ids=NaN(nc,4);
%for n=1:nc
%    tline=str2num(fgetl(fid));
%    ids(n,:)=tline(1:4);
%end
%Compute theta
theta=NaN(nc,1);
%n=19555;
for n=1:nc
    xs=xc(n);ys=yc(n);lats=latc(n);lons=lonc(n);
    x=NaN(3,1);y=NaN(3,1);lat=NaN(3,1);
    x(1)=xn(ids(n,2));x(2)=xn(ids(n,3));x(3)=xn(ids(n,4));
    y(1)=yn(ids(n,2));y(2)=yn(ids(n,3));y(3)=yn(ids(n,4));
    lat(1)=latn(ids(n,2));lat(2)=latn(ids(n,3));lat(3)=latn(ids(n,4));
    %lon(1)=lonn(ids(n,2));lon(2)=lonn(ids(n,3));lon(3)=lonn(ids(n,4));
    test=sign(lat-lats);
    i1=find(test==1);
    i2=find(test==-1);
    if length(i1)==2
        p1=i1;p2=i2;
    else
        p1=i2;p2=i1;
    end
    %Find positions at cell edge with lon=lons and x=xs
    xlat=NaN(2,1);
    ylat=NaN(2,1);
    xx=NaN(2,1);
    yx=ys*ones(2,1);
    %lons
    a=(lat(p1)-lat(p2))./(x(p1)-x(p2));
    b=lat(p1)-(lat(p1)-lat(p2))./(x(p1)-x(p2)).*x(p1);
    xlat=(lats-b)./a;
    a=(lat(p1)-lat(p2))./(y(p1)-y(p2));
    b=lat(p1)-(lat(p1)-lat(p2))./(y(p1)-y(p2)).*y(p1);
    ylat=(lats-b)./a;
    %y
    test=sign(y-ys);
    i1=find(test==1);
    i2=find(test==-1);
    if length(i1)==2
        p1=i1;p2=i2;
    else
        p1=i2;p2=i1;
    end
    a=(y(p1)-y(p2))./(x(p1)-x(p2));
    b=y(p1)-(y(p1)-y(p2))./(x(p1)-x(p2)).*x(p1);
    xx=(ys-b)./a;
    %Find sin(angle) by cross product.
    %ax,ay is the components of the vector defining the direction of the
    %x-axis, while bx,by gives the vector defining the lon-axis.
    %ay is by definition zero, which simplifies the expression
    ax=xx-xs;
    %the sign of the angle is given by -sign(by./bx)
    bx=xlat-xs;by=ylat-ys;
    theta_tmp=abs(asin(by./sqrt(bx.^2+by.^2).*sign(ax))).*(-sign(by./bx));
    theta(n)=mean(theta_tmp);
end
%Remove nans from theta. 
if any(isnan(theta))
  inan=find(isnan(theta));
  innan=find(~isnan(theta));
  for n=1:length(inan)
    ds2=(xc(innan)-xc(inan(n))).^2+(yc(innan)-yc(inan(n))).^2;
    imin=find(ds2==min(ds2));
    theta(inan(n))=theta(innan(imin));
  end
end

if 0
    p=plot(x,y,'k.');
    hold on
    p=plot(xs,ys,'ko');
    p=plot(xlat,ylat,'ro');
    p=plot(xx,yx,'rx');
end