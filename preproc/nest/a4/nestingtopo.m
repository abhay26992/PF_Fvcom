function Mobj = nestingtopo(hn,coef_n,ind_n,R,Mobj)

%%%%%%%%%nesting grid on the nudging zone%%%%%%%%%%%%%%
% load grd
% use grd

xn = Mobj.x;
yn = Mobj.y;
nn = length(xn);
disp('Interpolating h to node points on the nudging zone')

hnord = sum(coef_n.*hn(ind_n),2);



Rdist = NaN(nn,1);

[x_obc,y_obc,obcnodes] = get_obc(Mobj);
for n = 1:nn
    dist  = sqrt((x_obc-xn(n)).^2+(y_obc-yn(n)).^2);
    Rdist(n) = min(dist);
end

%Find weight function
r1 = 1.1*R;
r2 = 2*R;
weigth = NaN(Mobj.nVerts,1);
i = find(Rdist<r1);
weigth(i) = 1;
i = find(Rdist>r2);
weigth(i) = 0;
i =  find(Rdist >= r1 & Rdist <=r2);
a = 1/(r1-r2);
b = r2/(r2-r1);
weigth(i) = a.*Rdist(i)+b;
Mobj.hfv = Mobj.h;
hfv = Mobj.h;
h = hnord.*weigth+hfv.*(1-weigth);
Mobj.h = h;


