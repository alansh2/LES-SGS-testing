% stress.m ... 

clear
close all

% start parallel pool for Gaussian filter
pc = parcluster('local');
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));

% read data
load('data.mat')

[ny,nx,nz] = size(X);
x = X(1,:,1);
y = Y(1:ny);
z = reshape(Z(1,1,:),1,nz);

dx = x(2) - x(1);
dz = z(2) - z(1);

h = y(end);
L = x(end);
D = z(end);

% SGS parameters
P = 4*h + 2*D; % wetted perimeter
d_h = 8*h*D/P; % hydraulic diameter
% characteristic length scale as defined by
% https://www.cfd-online.com/Wiki/Hydraulic_diameter
Delta = 0.07*d_h;
Gstd = sqrt(Delta^2/12)*[1/dx 1 1/dz]; % standard deviation of Gaussian filter
save('properties.mat','Delta');

% resolved scales from Gaussian filter
Gu1 = GaussianFilter(y,U,Gstd);
Gu2 = GaussianFilter(y,V,Gstd);
Gu3 = GaussianFilter(y,W,Gstd);
save('resolved.mat','Gu1','Gu2','Gu3','-v7.3');
% SGS stress:   ___    _ _
%         Tij = uiuj = uiuj
tau = matfile('T.mat','Writable',true);
tau.T11 = GaussianFilter(y,U.*U,Gstd) - Gu1.*Gu1;
tau.T12 = GaussianFilter(y,U.*V,Gstd) - Gu1.*Gu2;
tau.T13 = GaussianFilter(y,U.*W,Gstd) - Gu1.*Gu3;
tau.T22 = GaussianFilter(y,V.*V,Gstd) - Gu2.*Gu2;
tau.T23 = GaussianFilter(y,V.*W,Gstd) - Gu2.*Gu3;
tau.T33 = GaussianFilter(y,W.*W,Gstd) - Gu3.*Gu3;
