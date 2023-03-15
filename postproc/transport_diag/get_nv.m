function nv=get_nv(nc)
%Usage: nv=get_nv(nc)
%Input: nc - number of cells.
%Output: nv(n,1:3) node id of three nodes in cell n. 
%disp('Reading grid file')

%casestr=get_cstr;
%fname_grd=['input/' casestr '_grd.dat'];

fid=fopen('C:/PhD/FVCOM/Matlab_Repository/Petermann_Bathy/Bathymetry_IceBridge/Input_Files/Natural_Sill/grd.dat','r');
fgetl(fid);
fgetl(fid);
nv=NaN*ones(nc,3);
for n=1:nc
  tline=str2num(fgetl(fid));
  cid=tline(2:4);
  nv(n,:)=cid;
end
