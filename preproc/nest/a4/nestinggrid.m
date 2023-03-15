function ngrd = nestinggrid(Rin,oend1,oend2,proj)
% Finds xc,yc,xn,yn,nv,lonn,latn,lonc and latc for the nestinggrid and puts it in the structure array ngrd. 
% Rin is the radius from the boundary or (if it is character) Rin is the filename of a
% nesting file created by fvcom. All cells within a radius of Rin from the boundary is within the nesting %zone
% oend1 and oend2 are logical variables. if oend1 is true, the distance that defines the nesting zone is calculated from the obc nodes which is a distance Rin from the start of the actial obc. if oend2 is true, the same is true at the end of the actual obc. The obc-nodes are always given in counter-clockwise direction.
% If proj is given ('utm' or 'pst') the angle to geographic east is also added to ngrd.
% Usage  ngrd=nestinggrid(Rin,proj)

createnew = 0;
if ~exist('./ngrd.mat','file')
   createnew = 1;
elseif ischar(Rin)
   load ngrd.mat
   if strcmp(Rin,ngrd.R)
      createnew = 1;
   end
else
   load ngrd.mat
   if ngrd.R ~= Rin;
      createnew = 1;
   end
end

if createnew
   grd = get_grd;
   use grd
   load M
   if ~ischar(Rin) 	% Find cells within a radius from the boundary
      [x_obc,y_obc,obcnodes] = get_obc(Mobj);
      if oend1
         x0 = x_obc(1);
	 y0 = y_obc(1);
         dist = sqrt((x_obc-x0).^2+(y_obc-y0).^2);
         i = find(dist >= Rin-0.1*Rin);
         x_obc = x_obc(i);
         y_obc = y_obc(i);
      end
      if oend2
         x0 = x_obc(end);
 	 y0 = y_obc(end);
         dist = sqrt((x_obc-x0).^2+(y_obc-y0).^2);
         i = find(dist >= Rin-0.1*Rin);
         x_obc = x_obc(i);
         y_obc = y_obc(i);
      end
      nestcells = [];
      for n = 1:nc
          dist = sqrt((x_obc-xc(n)).^2+(y_obc-yc(n)).^2);
          if min(dist) < Rin
             nestcells = [nestcells,n];
          end
      end
   end
   nv = get_nv(nc);
   if nargin == 4
      fvangle = anglec(Mobj);
      if any(isnan(fvangle))
         i1 = find(isnan(fvangle));
         i2 = find(~isnan(fvangle));
         F = TriScatteredInterp(xc(i2),yc(i2),fvangle(i2),'nearest');
         fvangle(i1) = F(xc(i1),yc(i1));
      end
   end
   %Nesting grid
   %Nodes
   if ~ischar(Rin)
      nid = [nv(nestcells,1),nv(nestcells,2),nv(nestcells,3)];
      nid = unique(nid);
      xn = xn(nid);
      yn = yn(nid);
   else
      ncid = netcdf.open(Rin,'nowrite');
      id = netcdf.inqVarID(ncid,'x');
      xnode = double(netcdf.getVar(ncid,id));
      id = netcdf.inqVarID(ncid,'y');
      ynode = double(netcdf.getVar(ncid,id));
      nid = [];
      grd = get_grd;
      use grd
      for n = 1:length(xnode)
          d = sqrt((yn-ynode(n)).^ 2+(xn-xnode(n)).^2);
          i = find(d == min(d));
          nid = [nid,i];
      end
      xn = xn(nid);
      yn = yn(nid);
      netcdf.close(ncid)
   end
   %h=get_depth;
   h = Mobj.h;
   h = h(nid);

   %Cells 
   % A second pass finding all the cells made from nodes in nid.
   %nestcells=[];
   %for n=1:nc                       
   %  if sum(ismember(nv(n,:),nid))==3
   %    nestcells=[nestcells;n];
   %  end
   %end
   xc = xc(nestcells);
   yc = yc(nestcells);

   nv = nv(nestcells,:);
   for n = 1:size(nv,1)
       for k = 1:3
           i = find(nid == nv(n,k));
           nv(n,k) = i;
       end
   end
   if nargin == 4
      fvangle = fvangle(nestcells);
   end

   %lat,lons
   %utmz='33 W';
   %utmzone=repmat(utmz,length(xn),1);
   %[latn,lonn]=utm2deg(xn,yn,utmzone);
   %utmzone=repmat(utmz,length(xc),1);
   %[latc,lonc]=utm2deg(xc,yc,utmzone);
   [lonn,latn] = my_project(xn,yn,'reverse');
   [lonc,latc] = my_project(xc,yc,'reverse'); 

   ngrd.xc = xc;
   ngrd.yc = yc;
   ngrd.xn = xn;
   ngrd.yn = yn;
   ngrd.nv = nv;
   ngrd.lonn = lonn;
   ngrd.latn = latn;
   ngrd.lonc = lonc;
   ngrd.latc = latc;
   ngrd.h = h; 
   if nargin == 4
      ngrd.fvangle = fvangle;
   end
   ngrd.R = Rin;
   ngrd.nid = nid;
   ngrd.cid = nestcells';%'
   ngrd.oend1 = oend1;
   ngrd.oend2 = oend2;
   save ngrd ngrd
   else
   load ngrd.mat
end
