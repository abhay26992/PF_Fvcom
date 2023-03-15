%VIndex  finds  and save the vertical interpolation index and coefficients for each fvcom node and cell
%center

clear all
close all




%%%%%%%%%%%%%ROMS grid%%%%%%%%%%%%%%%%
%saload('/global/work/qin/nclm_norkyst800/ROGrd.mat')
load  ROGrd
use  ROGrd

load M


load nearest4_r.mat
load nearest4_r2c.mat

disp('Vertical interpolation indices for fvcom cell and rho')
Vinterp=VInterpIndex(xr,yr,zr,Mobj,coef_r,ind_r,coef_r2c,ind_r2c);


save  -v7.3 Vinterp.mat  Vinterp


