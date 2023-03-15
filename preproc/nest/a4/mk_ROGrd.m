clear all
close all

%%% mk_ROGrd makes a structure of the ocean grid information from the roms. 

createnew=0;
if ~exist('./ROGrd.mat','file')
 createnew=1;
end

if createnew

%%%%%%%%Roms coordinates

datasource='arctic4_daily';

path='/global/work/apn/A4_modelruns/A4_nudging/linked_with_dates/';

createnew1=0;
if ~exist('./RCoor.mat','file')
 createnew1=1;
end

if createnew1
RCoor=RomsCoordinates(datasource,path);
save RCoor RCoor
else
load  RCoor 
end

%%%%%%%%%%%%find Roms Ocean grid%%%%%%%%
load rind

ROGrd=RomsOceanGrd(RCoor,rind);

save ROGrd ROGrd

else
    
    load ROGrd.mat
end
