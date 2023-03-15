function VIndex=VInterpIndex(xr,yr,zr,Mobj,coef_r,ind_r,coef_r2c,ind_r2c)


%VInterpIndex calculates coefficients and indices for fast vertical linear
%interpolation
% Input
%      xr,yr,zr: x-coordinate,y-coordinate and layer depth at roms rho ocean points
%      Mobj: fvcom mesh grid
%      coef_r,ind_r: coefficients and indices for horizontal interpolation
%      from roms rho point to fvcom node 
%      coef_r2c,ind_r2c: coefficients and indices for horizontal
%      interpolation from roms rho point to fvcom cell
% Output
%     VIndex has following fields
%      zcoef_n,zind_n: coefficients and indices for fvcom nodes
%      zcoef_c,zcoef_c: coefficients and indices for fvcom cell
%      ZN,ZC  : interpolated roms layer depths at fvcom node and cell
%      siglayz_c: fvcom sigma layer depth at cell center 
%
%

%%%first interpolate roms layer depth on roms grid to fvcom node and cell
%%%center

xn=Mobj.x;
yn=Mobj.y;
xc=Mobj.xc;
yc=Mobj.yc;
h=Mobj.h;
nn=length(xn);
nc=length(xc);
siglayz=Mobj.siglayz;
kb=get_kb;
siglayz_c=NaN(nc,kb-1);
nv=Mobj.tri;

for i=1:nc
  
  siglayz_c(i,:)=(siglayz(nv(i,1),:)+siglayz(nv(i,2),:)+siglayz(nv(i,3),:))/3;
  
end


%interpolating roms layer depths at rho points to nodes and cell center


[nx,nz]=size(zr);

ZN=NaN(nn,nz);
ZC=NaN(nc,nz);


for i=1:nz
    
   ZZ=zr(:,i);
   ZN(:,i)=sum(coef_r.*ZZ(ind_r),2); 
   ZC(:,i)=sum(coef_r2c.*ZZ(ind_r2c),2);
    
end


%%%looking for the index and coeffient for vertical interpolation at nodes
%%%and  cell centers

%node, for each point in siglayz, find its nearest layer depth points from
%roms, and corresponding index

zind_n=NaN(size(siglayz));
zind_n=repmat(zind_n,1,1,2);
zcoef_n=zind_n;

for i=1:nn 
    
    ZZ=ZN(i,:);
    
   for  j=1:kb-1
       
       SZ=siglayz(i,j);
       
       q1=find(ZZ>=SZ);
       q2=find(ZZ<SZ);
       
       if ~isempty(q1)&~isempty(q2)
        
           zind_n(i,j,:)=[q2(end) q1(1)];
           
           dz=[abs(ZZ(q2(end))-SZ) abs(ZZ(q1(1))-SZ)];
           dz=1./dz;
           
           zcoef_n(i,j,:)=dz./sum(dz);
           
       elseif isempty(q2)
           
           zind_n(i,j,:)=[q1(1) q1(1)];
           zcoef_n(i,j,:)=[0.5 0.5];
           
           
           
       else
           
           zind_n(i,j,:)=[q2(end) q2(end)];
           zcoef_n(i,j,:)=[0.5 0.5];
       end
       
    
   end     
end

%cell center

zind_c=NaN(size(siglayz_c));
zind_c=repmat(zind_c,1,1,2);
zcoef_c=zind_c;

for i=1:nc
    
    ZZ=ZC(i,:);
   
    
   for  j=1:kb-1
       
       SZ=siglayz_c(i,j);
       
       q1=find(ZZ>=SZ);
       q2=find(ZZ<SZ);
       
       if ~isempty(q1)&~isempty(q2)
        
           zind_c(i,j,:)=[q2(end) q1(1)];
           
           dz=[abs(ZZ(q2(end))-SZ) abs(ZZ(q1(1))-SZ)];
           dz=1./dz;
           
           zcoef_c(i,j,:)=dz./sum(dz);
           
       elseif isempty(q2)
           
           zind_c(i,j,:)=[q1(1) q1(1)];
           zcoef_c(i,j,:)=[0.5 0.5];
           
           
       else
           
           zind_c(i,j,:)=[q2(end) q2(end)];
           zcoef_c(i,j,:)=[0.5 0.5];
       end
   
   end     
end

make VIndex zcoef_n zind_n zcoef_c zind_c ZN ZC siglayz_c

 