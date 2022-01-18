%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A wavenumber-domain approach for quasi 3-D forward modeling of gravity anomalies and gradients.
% Author: Lianghui Guo (guolh@cugb.edu.cn), Yatong Cui
% Organization: China University of Geosciences (Beijing), School of Geophysics and Information Technology
% Compiled version: MATLAB R2017b
% Reference:
%       Cui Y T, Guo L H. A wavenumber-domain iterative approach for rapid 3-D imaging
%       of gravity anomalies and gradients. IEEE Access, 2019, 7: 34179-34188.
% Description of the input parameters: 
%       infile_msh: model mesh file
%       infile_mod: model density file, unit: g/cm^3
% Description of the output parameters: 
%       outfile_Uz: calculated anomaly
% Description of primary identifiers£º
%       x, y: x, y verctor
%       nx, ny: number of points in x and y directions
%       dx, dy: spacing in x and y directions
%       npts: extension points
%       p£ºdensity distribution, unit: g/cm^3
%       G£ºgravitational constant, unit: m^3/kg¡¤s^2
%       pex£ºdensity distribution after extension
%       Uz£ºgravity anomaly, unit: mGal
% Description of subroutine function: 
%       readmsh.m: read mesh file
%       readmod.m: read model file
%       wave2d.m: calculate wavenumber
%       extend_copy2d.m: copy edge extension
%       forward_Uz.m: calculated gravity anomaly
%       savegrd.m: save surfer text grd file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%%%%%%%%%%% I/O parameters %%%%%%%%%%%%%
infile_msh = 'model_stack5_101_101_41.msh'; 
infile_mod = 'model_stack5_101_101_41.mod';
outfile_Uz = 'Uz_2D.grd'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y,z,nx,ny,nz,dx,dy,dz] = readmsh(infile_msh);
p = readmod(infile_mod,nx,ny,nz);  p = p*1000;
nmax = max([nx ny]);
npts = 2^nextpow2(nmax);
pex = extend_copy2d(p,nx,ny,nz,npts);
G = 6.67e-11;  
% forward 
Uz = forward_Uz(pex,G,nx,ny,nz,dz,z,npts,dx,dy);
% save
savegrd(Uz,x,y,nx,ny,outfile_Uz);