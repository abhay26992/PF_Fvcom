% save_nearest4_r.m finds interpolation index and coefficients from ROMS
% rho points to fvcom nodes

load M
%%%%%%%%%%%%%ROMS Ocean grid%%%%%%%%%%%%%%%%

load ROGrd

xr=ROGrd.xr;
yr=ROGrd.yr;


disp('Roms rho point to fvcom node')
[coef_r,ind_r]=nearest4(Mobj,xr,yr,'node'); % from roms rho point to node

%save nearest4_r  coef_r ind_r
 save -v7.3  nearest4_r.mat coef_r ind_r
