function R=res_node(nc,xn,yn)
%Usage: R=res_node(nc,xn,yn)
%Input: nc - number of cells (
%       xn, yn - coordinates at node points
%Finds the resolution at each node as the distance to the nearest node.

%casestr=get_cstr;
casestr='pf_o_1';
fname_grd=['/home/abhay/Matlab_Repository/Petermann_Bathy/input/' casestr '_grd.dat'];
%Find shortest edge in each cell.
fid=fopen(fname_grd,'r');
nn=length(xn);
Rtmp=NaN*ones(nn,16);%One node can be included in up to 8 elements.
tline=fgetl(fid);
tline=fgetl(fid);
for n=1:nc
	if mod(n,100000)==0
	disp(n)
        end
  tline=str2num(fgetl(fid));
  nid=tline(2:4);
  d12=sqrt((xn(nid(1))-xn(nid(2))).^2+(yn(nid(1))-yn(nid(2))).^2);
  d13=sqrt((xn(nid(1))-xn(nid(3))).^2+(yn(nid(1))-yn(nid(3))).^2);
  d23=sqrt((xn(nid(3))-xn(nid(2))).^2+(yn(nid(3))-yn(nid(2))).^2);
  n1=sum(~isnan(squeeze(Rtmp(nid(1),:))));
  n2=sum(~isnan(squeeze(Rtmp(nid(2),:))));
  n3=sum(~isnan(squeeze(Rtmp(nid(3),:))));
  Rtmp(nid(1),n1+1:n1+2)=[d12,d13];
  Rtmp(nid(2),n2+1:n2+2)=[d12,d23];
  Rtmp(nid(3),n3+1:n3+2)=[d13,d23];
end
fclose(fid);
R=min(Rtmp,[],2);
