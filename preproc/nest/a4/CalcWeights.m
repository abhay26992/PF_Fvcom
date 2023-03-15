function [weight_cell,weight_node]=CalcWeights(time,w1,w2)
%Calculates linear weights in the nesting zone, from w1 at the obc to w2 
%at the inner end of the nesting zone. At the obc nodes, weights equals 1.
load M
[x_obc,y_obc,obcnodes]=get_obc(Mobj);
load ngrd
use ngrd
nc=length(xc);
nn=length(xn);
%Find max radius and node distance vector
if oend1
  x0=x_obc(1);y0=y_obc(1);
  dist=sqrt((x_obc-x0).^2+(y_obc-y0).^2);
  i=find(dist>R);
  x_obc=x_obc(i);
  y_obc=y_obc(i);
end  
if oend2
  x0=x_obc(end);y0=y_obc(end);
  dist=sqrt((x_obc-x0).^2+(y_obc-y0).^2);
  i=find(dist>R);
  x_obc=x_obc(i);
  y_obc=y_obc(i);
end
R=[];
d_node=NaN(nn,1);
for n=1:nn
    dist=min(sqrt((x_obc-xn(n)).^2+(y_obc-yn(n)).^2));
    d_node(n)=dist;
    R=[R,dist];
end
R=max(R);
%Interpolation values
DD=[0 R];
TT=[w1 w2];
%Cell distance vectors
d_cell=NaN(nc,1);
for n=1:nc
    dist=min(sqrt((x_obc-xc(n)).^2+(y_obc-yc(n)).^2));
    d_cell(n)=dist;
end
%Interpolate find weights
weight_node=interp1(DD,TT,d_node);
weight_cell=interp1(DD,TT,d_cell);
i=find(weight_node<0);
weight_node(i)=0;
%Use time to get weights in the right dimension for netcdf file
tdim=length(time);
weight_cell=repmat(weight_cell,1,tdim);
weight_node=repmat(weight_node,1,tdim);
