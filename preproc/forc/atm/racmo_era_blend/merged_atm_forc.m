function merged_out=merged_atm_forc(Mobj,lon,lat,RACMO_out,ERA_out,type)

% Write out merged atmospheric forcing using RACMO and ERA for the FVCOM
% grid. Usage:

% INPUT -> 
% 1) Mobj file of the latest FVCOM setup
% 2) longitude (lon) and latitude (lat) of the RACMO grid
% 3) Racmo output (RACMO_out) of the variable to be written out as .mat
% 4) Era output (ERA_out) of the same variable as .mat
% 5) Type of output : Defined over "node" or "cell"

% OUTPUT -> merged atmospheric input product of the variable for FVCOM 

%% Finding nodes and cells that fall outside the RACMO grid:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1a: Identify RACMO bounding box
rac_x_box=[lon(1:end,1); lon(end,1:end)'; lon(end:-1:1,end); lon(1,end:-1:1)'];
rac_y_box=[lat(1:end,1); lat(end,1:end)'; lat(end:-1:1,end); lat(1,end:-1:1)'];

% Step 2a: Find nodes in FVC domain outside RACMO

%nodes_outside_racmo = find(~inpolygon(Mobj.lon,Mobj.lat, rac_x_box, rac_y_box));
nodes_outside_racmo = ~inpolygon(Mobj.lon,Mobj.lat, rac_x_box, rac_y_box);
nodes_inside_racmo = inpolygon(Mobj.lon,Mobj.lat, rac_x_box, rac_y_box);

save rac_x_box rac_x_box
save rac_y_box rac_y_box
save nodes_outside_racmo nodes_outside_racmo
save nodes_inside_racmo nodes_inside_racmo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2b: 
[latc,lonc]=polarstereo_inv(Mobj.xc,Mobj.yc);
c_outside_racmo = ~inpolygon(lonc,latc, rac_x_box, rac_y_box);
c_inside_racmo = inpolygon(lonc,latc, rac_x_box, rac_y_box);

save c_outside_racmo c_outside_racmo
save c_inside_racmo c_inside_racmo
%% Merge the two products:

merged_out=NaN(length(RACMO_out),size(RACMO_out,2));
if type=="node"
    for i=1:length(RACMO_out)
        if nodes_inside_racmo(i,1)==1
             merged_out(i,:)=RACMO_out(i,:);
        else
            merged_out(i,:)=ERA_out(i,:);
        end
    end
elseif type=="cell"
   for i=1:length(RACMO_out)
       if c_inside_racmo(i,1)==1
           merged_out(i,:)=RACMO_out(i,:);
       else
           merged_out(i,:)=ERA_out(i,:);
       end
   end
end
           
    
