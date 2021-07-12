% testLES_SGS.m - LES SGS Model Test From DNS
% 
% DNS variables are roughly 11 GB in size.
% This script should use less than 100 GB at a time.

close all
clear all

% load DNS
dat = matfile('data.mat');

% specify grid widths
[ny, nx, nz] = size(dat,'U');
hx = 8*pi/(nx-1);
hy = 2/(ny-1);
hz = 3*pi/(nz-1);

delta = 1; % channel half height

% ---- Smagorinsky Model ----
Delta = (hx*hy*hz)^(1/3);
cs = 0.16;

% Gaussian filter
Gstd = delta/2; % Delta/2;
Gu1 = imgaussfilt3(dat.U,Gstd);
Gu2 = imgaussfilt3(dat.V,Gstd);
Gu3 = imgaussfilt3(dat.W,Gstd);

% calculate gradients of filtered data
[dUdy,dUdx,dUdz] = gradient(Gu1,hy,hx,hz);
save('gradient.mat','dUdx','dUdy','dUdz','-v7.3');
clear dUdx dUdy dUdz
grad = matfile('gradient.mat','Writable',true);
[dVdy,dVdx,dVdz] = gradient(Gu2,hy,hx,hz);
grad.dVdx = dVdx;
grad.dVdy = dVdy;
grad.dVdz = dVdz;
clear dVdx dVdy dVdz
[dWdy,dWdx,dWdz] = gradient(Gu3,hy,hx,hz);
grad.dWdx = dWdx;
grad.dWdy = dWdy;
grad.dWdz = dWdz;
clear dWdx dWdy dWdz

% SGS stress tensor = <uiuj> - <ui><uj>
% <.> = filtering operation
T11 = imgaussfilt3(dat.U.*dat.U,Gstd) - Gu1.*Gu1;
save('SGSstress.mat','T11','-v7.3');
clear T11
SGS = matfile('SGSstress.mat','Writable',true);
SGS.T12 = imgaussfilt3(dat.U.*dat.V,Gstd) - Gu1.*Gu2;
SGS.T13 = imgaussfilt3(dat.U.*dat.W,Gstd) - Gu1.*Gu3;
SGS.T22 = imgaussfilt3(dat.V.*dat.V,Gstd) - Gu2.*Gu2;
SGS.T23 = imgaussfilt3(dat.V.*dat.W,Gstd) - Gu2.*Gu3;
SGS.T33 = imgaussfilt3(dat.W.*dat.W,Gstd) - Gu3.*Gu3; % symmetric
clear Gu1 Gu2 Gu3

% strain rate tensor
S11 = grad.dUdx;
S12 = (grad.dUdy + grad.dVdx)/2;
S13 = (grad.dUdz + grad.dWdx)/2;
S22 = grad.dVdy;
S23 = (grad.dVdz + grad.dWdy)/2;
S33 = grad.dWdz; % symmetric

% calculate invariant
Snorm = sqrt(2*(S11.*S11 + 2*S12.*S12 + 2*S13.*S13 + S22.*S22 + 2*S23.*S23 + S33.*S33));

% stress tensor
nu = -2*(cs*Delta)^2*Snorm; % eddy viscosity
clear Snorm
T11 = nu.*S11;
save('SmagStress.mat','T11','-v7.3');
clear T11
smag = matfile('SmagStress.mat','Writable',true);
smag.T12 = nu.*S12;
smag.T13 = nu.*S13;
smag.T22 = nu.*S22;
smag.T23 = nu.*S23;
smag.T33 = nu.*S33;
clear S11 S12 S13 S22 S23 S33

% % if we want to load into variables
% x1 = dat.X;
% x2 = dat.Y;
% x3 = dat.Z;

% % ---- Plotting ----
% figure(1), clf
% contourf(dat.X(1:ny,1:nx,1), dat.Y(1:ny,1:nx,1), SGS.T12(1:ny,1:nx,1), 30, 'LineStyle', 'none');
% colorbar
% 
% figure(2), clf
% contourf(dat.X(1:ny,1:nx,1), dat.Y(1:ny,1:nx,1), smag.T12(1:ny,1:nx,1), 30, 'LineStyle', 'none');
% colorbar;

% % ---- A Priori ----
% AB11 = T11.*Tsmag11 + T12.*Tsmag12 + T13.*Tsmag13;
% AB22 = T12.*Tsmag12 + T22.*Tsmag22 + T23.*Tsmag23;
% AB33 = T13.*Tsmag13 + T23.*Tsmag23 + T33.*Tsmag33;
% AA11 = T11.*T11 + T12.*T12 + T13.*T13;
% AA22 = T12.*T12 + T22.*T22 + T23.*T23;
% AA33 = T13.*T13 + T23.*T23 + T33.*T33;
% BB11 = Tsmag11.*Tsmag11 + Tsmag12.*Tsmag12 + Tsmag13.*Tsmag13;
% BB22 = Tsmag12.*Tsmag12 + Tsmag22.*Tsmag22 + Tsmag23.*Tsmag23;
% BB33 = Tsmag13.*Tsmag13 + Tsmag23.*Tsmag23 + Tsmag33.*Tsmag33;
% trAB = AB11 + AB22 + AB33;
% trAA = AA11 + AA22 + AA33;
% trBB = BB11 + BB22 + BB33;
% angles = acosd(trAB./sqrt(trAA.*trBB));
