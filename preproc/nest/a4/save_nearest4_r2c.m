% save_nearest4_r.m finds interpolation index and coefficients from ROMS
% rho points to fvcom cell

load M
%%%%%%%%%%%%%ROMS grid%%%%%%%%%%%%%%%%
load  ROGrd

xr=ROGrd.xr;
yr=ROGrd.yr;

disp('Roms rho point to fvcom cell')
[coef_r2c,ind_r2c]=nearest4(Mobj,xr,yr,'cell'); % from roms rho point to cell center

save -v7.3 nearest4_r2c.mat coef_r2c ind_r2c
