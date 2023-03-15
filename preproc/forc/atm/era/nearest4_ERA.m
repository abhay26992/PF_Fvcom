clear all;close all;clc;

load ./eragrd
load ./M
load ./erawdind % index for the local domain
use erawdind

%xg and yg: RACMO x,y (converted) coordinates 
xg=xg(idr,jdr); %crop using idr,jdr
yg=yg(idr,jdr);

xr=xg(:);
yr=yg(:);



%%%%set those indices and coeffients as NaN if it is out of 2km wind domain
disp('indicies and coefficients for fvcom node')
xf=Mobj.x;
yf=Mobj.y;

nn=length(xf);
index=NaN(nn,4);  

 coef=NaN(nn,4);
      
for n=1:nn
   
   
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
     
N4ERA.ncoef=coef;
N4ERA.nindex=index;

clear coef index xf yf n nn IN


%%%%%%%%%%%%%%%%%indicies and coefficients for cell center%%%%%%%

disp('indicies and coefficients for fvcom cell')
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
    %%%if fvcom cell is surrounded by four points, then the index are the
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
       
       N4ERA.ccoef=coef;
       N4ERA.cindex=index;
       
   
    
save -v7.3  N4ERA N4ERA
