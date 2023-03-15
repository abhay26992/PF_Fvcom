function R=res_cell(nc,xn,yn)
%Usage: R=res_cell(nc,xn,yn)
%Input: nc - number of cells.
% xn, yn - x and y coordinates at node points
%res_cell finds the resolution, R, given by the length of the shortest edge in each cell.
disp('Reading grid file')

%[nc,nn]=GridSize;
%[xc,yc,xn,yn]=readgrid(nc,nn);

%casestr=get_cstr;
%fname_grd=['input/' casestr '_grd.dat'];
%Find shortest edge in each cell.
fid=fopen('/home/abhay/Matlab_Repository/Petermann_Bathy/input/pf_o_1_grd.dat','r');
disp('Find length of shortest edge in each cell')
R=NaN*ones(nc,1);
tline=fgetl(fid);
tline=fgetl(fid);
for n=1:nc
  tline=str2num(fgetl(fid));
  cid=tline(2:4);
  d1=sqrt((xn(cid(1))-xn(cid(2))).^2+(yn(cid(1))-yn(cid(2))).^2);
  d2=sqrt((xn(cid(1))-xn(cid(3))).^2+(yn(cid(1))-yn(cid(3))).^2);
  d3=sqrt((xn(cid(3))-xn(cid(2))).^2+(yn(cid(3))-yn(cid(2))).^2);
  R(n)=min([d1,d2,d3]);
end
