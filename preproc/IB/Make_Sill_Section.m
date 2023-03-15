% Make section plots to visualise the sill 

% Run after running Gen_Sill_Geom.m

%% Along the 3 lines:

%%%%%%%%%%%%%%%%%%%%%%%%% Center Line: Line 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h2=zeros(length(ib_lat2),1);
h2_old=zeros(length(ib_lat2),1);
cs_ind=zeros(length(ib_lat2),1);
h2_lat=zeros(length(ib_lat2),1);
h2_lon=zeros(length(ib_lat2),1);

for i=1:length(ib_lat2)
    l=spheredist(Mobj.lat,Mobj.lon,ib_lat2(i,1),ib_lon2(i,1));
    cs_ind(i,1)=find(l==min(l));
    h2(i,1)=new_h(cs_ind(i,1));
    h2_old(i,1)=Mobj.h(cs_ind(i,1));
    h2_lat(i,1)=Mobj.lat(cs_ind(i,1));
    h2_lon(i,1)=Mobj.lon(cs_ind(i,1));
end


h2_ind =find(inpolygon(h2_lon,h2_lat,xp,yp));
h2_lat=h2_lat(h2_ind);h2_lon=h2_lon(h2_ind);h2=h2(h2_ind);h2_old=h2_old(h2_ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Left: Line 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1=zeros(length(ib_lat1),1);
h1_old=zeros(length(ib_lat1),1);
cs_ind=zeros(length(ib_lat1),1);
h1_lat=zeros(length(ib_lat1),1);
h1_lon=zeros(length(ib_lat1),1);

for i=1:length(ib_lat1)
    l=spheredist(Mobj.lat,Mobj.lon,ib_lat1(i,1),ib_lon1(i,1));
    cs_ind(i,1)=find(l==min(l));
    h1(i,1)=new_h(cs_ind(i,1));
    h1_old(i,1)=Mobj.h(cs_ind(i,1));
    h1_lat(i,1)=Mobj.lat(cs_ind(i,1));
    h1_lon(i,1)=Mobj.lon(cs_ind(i,1));
end

h1_ind =find(inpolygon(h1_lon,h1_lat,xp,yp));
h1_lat=h1_lat(h1_ind);h1_lon=h1_lon(h1_ind);h1=h1(h1_ind);h1_old=h1_old(h1_ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right: Line 2f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h2_f=zeros(length(ib_lat2_f),1);
h2_f_old=zeros(length(ib_lat2_f),1);
cs_ind=zeros(length(ib_lat2_f),1);
h2_f_lat=zeros(length(ib_lat2_f),1);
h2_f_lon=zeros(length(ib_lat2_f),1);

for i=1:length(ib_lat2_f)
    l=spheredist(Mobj.lat,Mobj.lon,ib_lat2_f(i,1),ib_lon2_f(i,1));
    cs_ind(i,1)=find(l==min(l));
    h2_f(i,1)=new_h(cs_ind(i,1));
    h2_f_old(i,1)=Mobj.h(cs_ind(i,1));
    h2_f_lat(i,1)=Mobj.lat(cs_ind(i,1));
    h2_f_lon(i,1)=Mobj.lon(cs_ind(i,1));
end

h2_f_ind =find(inpolygon(h2_f_lon,h2_f_lat,xp,yp));
h2_f_lat=h2_f_lat(h2_f_ind);h2_f_lon=h2_f_lon(h2_f_ind);h2_f=h2_f(h2_f_ind);h2_f_old=h2_f_old(h2_f_ind);

%% Plot:

figure (1);clf      %new
scatter(h2_lon,h2_lat,100,h2);colorbar;colormap jet;caxis([550 850]);hold on;
scatter(h1_lon,h1_lat,100,h1);colorbar;colormap jet;
scatter(h2_f_lon,h2_f_lat,100,h2_f);colorbar;colormap jet;
scatter(Mobj.lon,Mobj.lat,0.0005,'k');

figure (2);clf      %old
scatter(h2_lon,h2_lat,100,h2_old);colorbar;colormap jet;caxis([550 850]);hold on;
scatter(h1_lon,h1_lat,100,h1_old);colorbar;colormap jet;
scatter(h2_f_lon,h2_f_lat,100,h2_f_old);colorbar;colormap jet;
scatter(Mobj.lon,Mobj.lat,0.0005,'k');

%% Cross - Section:

% Using figure (1) or (2):

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 298.56;
y1 = 80.83;
x2 = 299.4;
y2 = 80.9;
xn = 450;
xx = linspace(x1,x2,xn);
yy = linspace(y1,y2,xn);

%clear is
for i = 1:numel(xx)
    dum = spheredist(yy(i),xx(i)-360,Mobj.lat,Mobj.lon-360);
    is(i) = find(dum == min(dum));
end
is = unique(is,'stable');

h_cross1=new_h(is);h_crosslat_1=Mobj.lat(is);h_crosslon_1=Mobj.lon(is);
h_cross1_old=Mobj.h(is);

% For cross-section
zcrs = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-new_h),Mobj.siglay),Mobj.zisf); %change to old_h : use Mobj.h instead of new_h
zcrs = zcrs(is,:);
ll = repmat(h_crosslon_1,[1 23]);
dval=ones(size(ll));

figure(100);clf
pcolorjw(ll,-zcrs,dval);
set(gca,'Color','k');
ylim([-1000 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 298.76;
y1 = 80.73;
x2 = 299.59;
y2 = 80.83;
xn = 450;
xx = linspace(x1,x2,xn);
yy = linspace(y1,y2,xn);

%clear is2
for i = 1:numel(xx)
    dum = spheredist(yy(i),xx(i)-360,Mobj.lat,Mobj.lon-360);
    is2(i) = find(dum == min(dum));
end
is2 = unique(is2,'stable');

h_cross2=new_h(is2);h_crosslat_2=Mobj.lat(is2);h_crosslon_2=Mobj.lon(is2);
h_cross2_old=Mobj.h(is2);

% For cross-section
zcrs = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-new_h),Mobj.siglay),Mobj.zisf); %change to old_h : use Mobj.h instead of new_h
zcrs = zcrs(is2,:);
ll = repmat(h_crosslon_2,[1 23]);
dval=ones(size(ll));

figure(100);clf
pcolorjw(ll,-zcrs,dval);
set(gca,'Color','k');
ylim([-1000 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 299;
y1 = 80.66;
x2 = 299.72;
y2 = 80.77;
xn = 450;
xx = linspace(x1,x2,xn);
yy = linspace(y1,y2,xn);

for i = 1:numel(xx)
    dum = spheredist(yy(i),xx(i)-360,Mobj.lat,Mobj.lon-360);
    is3(i) = find(dum == min(dum));
end
is3 = unique(is3,'stable');

h_cross3=new_h(is3);h_crosslat_3=Mobj.lat(is3);h_crosslon_3=Mobj.lon(is3);
h_cross3_old=Mobj.h(is3);

% For cross-section
zcrs = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-new_h),Mobj.siglay),Mobj.zisf); %change to old_h : use Mobj.h instead of new_h
zcrs = zcrs(is3,:);
ll = repmat(h_crosslon_3,[1 23]);
dval=ones(size(ll));

figure(100);clf
pcolorjw(ll,-zcrs,dval);
set(gca,'Color','k');
ylim([-1000 0]);

%% Plot

figure (3);clf        %new
scatter(h_crosslon_1,h_crosslat_1,100,h_cross1);colorbar;colormap jet;caxis([550 850]);hold on;
scatter(h_crosslon_2,h_crosslat_2,100,h_cross2);colorbar;colormap jet;
scatter(h_crosslon_3,h_crosslat_3,100,h_cross3);colorbar;colormap jet;
scatter(Mobj.lon,Mobj.lat,0.0005,'k');

figure (4);clf        %old
scatter(h_crosslon_1,h_crosslat_1,100,h_cross1_old);colorbar;colormap jet;caxis([550 850]);hold on;
scatter(h_crosslon_2,h_crosslat_2,100,h_cross2_old);colorbar;colormap jet;
scatter(h_crosslon_3,h_crosslat_3,100,h_cross3_old);colorbar;colormap jet;
scatter(Mobj.lon,Mobj.lat,0.0005,'k');



%% All plots:

% new
figure (5);clf  
scatter(h2_lon,h2_lat,100,h2);colorbar;colormap jet;hold on;
scatter(h1_lon,h1_lat,100,h1);colorbar;colormap jet;
scatter(h2_f_lon,h2_f_lat,100,h2_f);colorbar;colormap jet;
scatter(h_crosslon_1,h_crosslat_1,100,h_cross1);colorbar;colormap jet;
scatter(h_crosslon_2,h_crosslat_2,100,h_cross2);colorbar;colormap jet;
scatter(h_crosslon_3,h_crosslat_3,100,h_cross3);colorbar;colormap jet;
caxis([550 850]);scatter(Mobj.lon,Mobj.lat,0.0005,'k');


%old
figure (6);clf      
scatter(h2_lon,h2_lat,100,h2_old);colorbar;colormap jet;hold on;
scatter(h1_lon,h1_lat,100,h1_old);colorbar;colormap jet;
scatter(h2_f_lon,h2_f_lat,100,h2_f_old);colorbar;colormap jet;
scatter(h_crosslon_1,h_crosslat_1,100,h_cross1_old);colorbar;colormap jet;
scatter(h_crosslon_2,h_crosslat_2,100,h_cross2_old);colorbar;colormap jet;
scatter(h_crosslon_3,h_crosslat_3,100,h_cross3_old);colorbar;colormap jet;
caxis([550 850]);scatter(Mobj.lon,Mobj.lat,0.0005,'k');

%% For cross - section like along section plots:

% Using figure (1) or (2):

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 298.56;
y1 = 80.89;
x2 = 299.23;
y2 = 80.65;
xn = 450;
xx = linspace(x1,x2,xn);
yy = linspace(y1,y2,xn);

%clear is
for i = 1:numel(xx)
    dum = spheredist(yy(i),xx(i)-360,Mobj.lat,Mobj.lon-360);
    is(i) = find(dum == min(dum));
end
is = unique(is,'stable');

h_along_1=new_h(is);h_alonglat_1=Mobj.lat(is);h_alonglon_1=Mobj.lon(is);
h_along1_old=Mobj.h(is);

% For cross-section like along-section plots
zcrs = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-Mobj.h),Mobj.siglay),Mobj.zisf); %change to old_h : use Mobj.h instead of new_h
zcrs = zcrs(is,:);
ll = repmat(h_alonglon_1,[1 23]);
dval=ones(size(ll));

figure(100);clf
pcolorjw(ll,-zcrs,dval);
set(gca,'Color','k');
ylim([-1000 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 298.85;
y1 = 80.93;
x2 = 299.56;
y2 = 80.68;
xn = 450;
xx = linspace(x1,x2,xn);
yy = linspace(y1,y2,xn);

%clear is2
for i = 1:numel(xx)
    dum = spheredist(yy(i),xx(i)-360,Mobj.lat,Mobj.lon-360);
    is2(i) = find(dum == min(dum));
end
is2 = unique(is2,'stable');

h_along_2=new_h(is2);h_alonglat_2=Mobj.lat(is2);h_alonglon_2=Mobj.lon(is2);
h_along2_old=Mobj.h(is2);

% For cross-section
zcrs = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-new_h),Mobj.siglay),Mobj.zisf); %change to old_h : use Mobj.h instead of new_h
zcrs = zcrs(is2,:);
ll = repmat(h_alonglon_2,[1 23]);


dval=ones(size(ll));

figure(100);clf
pcolorjw(ll,-zcrs,dval);
set(gca,'Color','k');
ylim([-1000 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line 3 (i.e. 2f) %%%%%%%%%%%%%%%%%%%%%%%%

x1 = 299.07;
y1 = 80.98;
x2 = 299.7;
y2 = 80.73;
xn = 450;
xx = linspace(x1,x2,xn);
yy = linspace(y1,y2,xn);

for i = 1:numel(xx)
    dum = spheredist(yy(i),xx(i)-360,Mobj.lat,Mobj.lon-360);
    is3(i) = find(dum == min(dum));
end
is3 = unique(is3,'stable');

h_along_3=new_h(is3);h_alonglat_3=Mobj.lat(is3);h_alonglon_3=Mobj.lon(is3);
h_along3_old=Mobj.h(is3);

% For cross-section
zcrs = bsxfun(@plus,bsxfun(@times,(Mobj.zisf-new_h),Mobj.siglay),Mobj.zisf); %change to old_h : use Mobj.h instead of new_h
zcrs = zcrs(is3,:);
ll = repmat(h_alonglon_3,[1 23]);
dval=ones(size(ll));

figure(100);clf
pcolorjw(ll,-zcrs,dval);
set(gca,'Color','k');
ylim([-1100 0]);