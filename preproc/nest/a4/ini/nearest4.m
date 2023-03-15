function  [coef, index]=nearest4(Mobj,xr,yr,invar)

% for each fvcom node or cell center, find its four nearest points from retangle grid 
% in order to perform fast linear interpolation
% INPUT
%    Mobj: Matlab mesh object
%    xr: x-coordinate of the retangle grid
%    yr: y-corrdinate of the retangle grid
%    invar: two options:{'node','cell'}
% OUTPUT
%    coef: a nx4 matrix of four coeffients for each node or cell point
%    index: a nx4 matrix of four index of the latitude and longitude
% NOTE
%    In this routine, we only use ocean points from the retangle grid. If
%    there are less than four points around it, indicating the point closes
%    to land, we choose the nearest grid as its representative. This means
%    that the index will be the same for four points, each of coefficient
%    of 0.25. 


switch invar
    

    case{'node'}
        
        xf=Mobj.x;
        yf=Mobj.y;
        nn=length(xf);
        %%%%find 4 nearest roms points to a fvcom point

       index=NaN(nn,4);

       coef=NaN(nn,4);

       for n=1:nn;
    
   
        dx=xr-xf(n);
        dy=yr-yf(n);
       ds=dx.^2+dy.^2;
    
       [ds0,ii]=sort(ds);
    
        tind=ii(1:4);
    
        xx=xr(tind)-repmat(xf(n),4,1);
    
        yy=yr(tind)-repmat(yf(n),4,1);
    
        clear i
    
        an=angle(xx+i*yy);
     
        [an1,ZZ]=sort(an);
    
        tind=tind(ZZ);
    
        xx=xr(tind)-repmat(xf(n),4,1);
    
        yy=yr(tind)-repmat(yf(n),4,1);
    
        xx(end+1)=xx(1);
        yy(end+1)=yy(1);
    
        for k=1:4
         v1=[xx(k);yy(k)];
         v2=[xx(k+1);yy(k+1)];     
         AN(k)=atan2(abs(det([v1,v2])),dot(v1,v2));     
         end
    
    %%%%%%%%%%%%%%%%%%%
    %%%if fvcom node is surrounded by four points, then the index are the
    %%%four points; otherwise, it indicates the node is close to land or
    %%%island, we choose the nearest point as index.
    
         if  abs(sum(AN)-2*pi)<1e-6 
            index(n,:)=tind;
            dx=xr(index(n,:))-xf(n);
            dy=yr(index(n,:))-yf(n);
            ds=sqrt(dx.^2+dy.^2);
            ds=1./ds;
            coef(n,:)=ds(:)./sum(ds);
         else
              index(n,:)=ii(1);
              coef(n,:)=1/4;
         end
    
     end
        
        
    case{'cell'}
       xf=Mobj.xc;
       yf=Mobj.yc;
       nn=length(xf);
        index=NaN(nn,4);

       coef=NaN(nn,4);

       for n=1:nn;
    
   
        dx=xr-xf(n);
        dy=yr-yf(n);
       ds=dx.^2+dy.^2;
    
       [ds0,ii]=sort(ds);
    
        tind=ii(1:4);
    
        xx=xr(tind)-repmat(xf(n),4,1);
    
        yy=yr(tind)-repmat(yf(n),4,1);
    
        clear i
    
        an=angle(xx+i*yy);
     
        [an1,ZZ]=sort(an);
    
        tind=tind(ZZ);
    
        xx=xr(tind)-repmat(xf(n),4,1);
    
        yy=yr(tind)-repmat(yf(n),4,1);
    
        xx(end+1)=xx(1);
        yy(end+1)=yy(1);
    
        for k=1:4
         v1=[xx(k);yy(k)];
         v2=[xx(k+1);yy(k+1)];     
         AN(k)=atan2(abs(det([v1,v2])),dot(v1,v2));     
         end
    
    %%%%%%%%%%%%%%%%%%%
    %%%if fvcom node is surrounded by four points, then the index are the
    %%%four points; otherwise, it indicates that the node is close to land or
    %%%island, we choose the nearest point as index, but set the cofficient to zero.
    
         if  abs(sum(AN)-2*pi)<1e-6 
            index(n,:)=tind;
            dx=xr(index(n,:))-xf(n);
            dy=yr(index(n,:))-yf(n);
            ds=sqrt(dx.^2+dy.^2);
            ds=1./ds;
            coef(n,:)=ds(:)./sum(ds);
         else
              index(n,:)=ii(1);
             % coef(n,:)=1/4;
             coef(n,:)=0;
         end
       end
       
end