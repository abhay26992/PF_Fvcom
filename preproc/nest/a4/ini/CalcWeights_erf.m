function [weight_cell,weight_node]=CalcWeights_erf(time,R)
%Calculates nonlinear weights in the nesting zone, from 1 at the obc to
% 0 at the end of nesting zone.  At the obc nodes, weights equals 1



load M
[x_obc,y_obc,obcnodes]=get_obc(Mobj);
load ngrd
use ngrd
nc=length(xc);
nn=length(xn);
%Find max radius
d_cell=NaN(nc,1);
for n=1:nc
    dist=min(sqrt((x_obc-xc(n)).^2+(y_obc-yc(n)).^2));
    d_cell(n)=dist;
    
end

xmid=R/2;
xwdth=R/7;
weight_cell=1-0.5.*(1+erf((d_cell-xmid)./(xwdth*sqrt(2))));

i=find(weight_cell<0);
weight_cell(i)=0;

d_node=NaN(nn,1);
for n=1:nn
    dist=min(sqrt((x_obc-xn(n)).^2+(y_obc-yn(n)).^2));
    d_node(n)=dist;
end

weight_node=1-0.5.*(1+erf((d_node-xmid)./(xwdth*sqrt(2))));

i=find(weight_node<0);
weight_node(i)=0;

%%%%%%%%%%%%%%%%set weights on obc nodes as 1%%%%%%%%%%%%%%%%%%%%%%

load tge




q1=find(isonb==2); % nodes on the open boundary

q2=find(isbce==2); % cells on the open boundary

weight_node(nid(q1),:)=1;
weight_cell(cid(q2),:)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




tdim=length(time);
weight_cell=repmat(weight_cell,1,tdim);
weight_node=repmat(weight_node,1,tdim);
