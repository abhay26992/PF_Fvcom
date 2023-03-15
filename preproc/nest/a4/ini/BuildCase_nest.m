clear all; close all;

global ftbverbose
ftbverbose = true; %print information to screen [true/false]

load M

%%%%%%%%%%%%%making nesting grid%%%%%%%%%%%%%%%%%%%%%%

load nearest4_r

load ROGrd

%casename=get_cstr;

casename='pf_o_1'

R=5000;



% interpolating roms topography


Mobj=nestingtopo(ROGrd.hr,coef_r,ind_r,R,Mobj);



save M Mobj

%Create new nesting-grid file with correct h
%!rm ngrd.mat
%ngrd=nestinggrid(nestfile,'pst');
ngrd=nestinggrid(R,1,1,'utm');

save ngrd ngrd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Mobj=read_sigma(Mobj,['/home/abhay/Matlab_Repository/Petermann_Bathy/input/' casename '_sigma.dat']);

save M Mobj


% dump mesh and connectivity
%write_FVCOM_grid(Mobj,['input/' casename '_grd.dat'])
% dump bathymetry
write_FVCOM_bath(Mobj,['/home/abhay/Matlab_Repository/Petermann_Bathy/input/' casename '_dep.dat'])
% % dump open boundary node list
% write_FVCOM_obc(Mobj,['input/' casename '_obc.dat'])
% dump sponge layer file
%write_FVCOM_sponge(Mobj,['input/' casename '_spg.dat']);
%% dump Coriolis file
%write_FVCOM_cor(Mobj,['input/' casename '_cor.dat'])



save M Mobj
