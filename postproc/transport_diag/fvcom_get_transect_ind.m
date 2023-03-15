function ind = fvcom_get_transect_ind(Mobj,xline,yline,dscale,type)
% Get the indices in a roms grid along the line between points in xl,
% lat_lim
tge
if nargin < 5
  type='cell';
end
if nargin <4
  dscale=50000;
end
if strcmp(type,'cell')
  x=Mobj.xc;y=Mobj.yc;
elseif strcmp(type,'node')
  x=Mobj.x;y=Mobj.y;
else
  disp(['Illegal type : ' type])
  return
end
% Model domain limits
x_min = min(x); x_max = max(x);
y_min = min(y); y_max = max(y);

% Check if requested transect is inside model domain
outside=0;
if xline(1)<x_min
    outside=1;
end
if xline(2)>x_max
    outside=1;
end
if yline(1)<y_min
    outside=1;
end
if yline(2)>y_max
    outside=1;
end

if outside
    disp(['Part of the requested transect lies outside the model domain. Please select new points!'])
    return
end



% Work in utm-coordinates (maybe use different projections for other cases)
%m_proj('utm', 'lon', [lon_min lon_max], 'lat', [lat_min lat_max], 'zon', 33, 'e%ll', 'wgs84');

% Plot map of grid area
%m_gshhs_f('patch',[.7 .7 .7]);
%m_grid


% Coordinates of grid points
%[x,y] = m_ll2xy(lon, lat);
%[x,y]=deg2utm(lat,lon,'33 W');

% Coordinates of transect endpoints
%[c1, c2] = m_ll2xy(lon_lim, lat_lim);
%[c1, c2] = deg2utm(lat_lim, lon_lim,'33 W');

distance1 = ((x - xline(1)).^2) + ((y - yline(1)).^2); % Distance from first point to all gridpoints 
ind1 = find(distance1 == min(min(distance1))); % Find indices of the first point
x1 = x(ind1); y1 = y(ind1); % Nearest grid grid point of first point

distance2 = ((x - xline(2)).^2) + ((y - yline(2)).^2); % Distance from second point to all gridpoints
ind2 = find(distance2 == min(min(distance2))); % Find indices of the second point
x2 = x(ind2); y2 = y(ind2); % Nearest grid grid point of second point

if x1>x2 
    x1_temp = x1;
    x1 = x2;
    x2 = x1_temp;
    y1_temp = y1;
    y1 = y2;
    y2 =y1_temp;
    ind1_temp = ind1;
    ind1 = ind2;
    ind2 = ind1_temp;
    clear x1_temp y1_temp;
end
 
% Find the equation for the straight line between the points 
a = (y2 - y1) / (x2 - x1); % slope
b = y1 - a * x1; % intersection

% Find indices along the line
ind=ind1;
while 1
 tmpind=nbe(ind(end),:);
 i=find(tmpind==0);
 tmpind(i)=[];
 xn=x(tmpind);
 yn=y(tmpind);
if length(ind)==1
  i=find((xn<x1)|(ismember(tmpind',ind)));%'
else
  i=find((ismember(tmpind',ind)));%'
end
 %xn(i)=[];yn(i)=[];
 ds2l=mean(abs(a.*xn+b-yn));
 % s=(xn-x2).^2+(yn-y2).^2;
 s1=abs(a.*xn+b-yn)-mean(abs(a.*xn+b-yn));
 s2=sqrt((xn-x2).^2+(yn-y2).^2)-mean(sqrt((xn-x2).^2+(yn-y2).^2));
 s=ds2l/dscale*s1+s2;
 if ~isempty(i);
  s(i)=[];
  tmpind(i)=[];
 end
 i=find(s==min(s));
 ind=[ind,tmpind(i)];
 if ind(end)==ind2
   break
 end
end
ind=ind';
%s = sqrt((x2 - x1)^2 + (y2 - y1)^2);
%theta = asin((y2 - y1) / s);
%
%dx = ds * cos(theta);
%xvar = x1:dx:x2;
%yvar = a * xvar + b;
%
%% Find indices along the line 
%
%k=1;
%ind = []; 
%for i=1:length(xvar)
%    distance = ((x-xvar(i)).^2) + ((y-yvar(i)).^2); % Distance from line to all gridpoints 
%    ii = find(distance == min(min(distance)));
%    if isempty(intersect(ind, ii))
%        ind(k) = ii;
%        k = k + 1;
%    end
%end

% Plot selected positions on map to check that they are correct
%hold on
%for n=1:length(indI)
%    m_plot(lon_rho(indI(n),indJ(n)),lat_rho(indI(n),indJ(n)),'.');
%end






