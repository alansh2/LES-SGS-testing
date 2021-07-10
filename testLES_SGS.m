% testLES_SGS.m - LES SGS Model Test From DNS
% 
% 

close all
clear all

% load DNS
dat = matfile('data.mat');

u1 = dat.U;
u2 = dat.V;
u3 = dat.W;

% Gaussian filter
Gstd = 0.5; % half channel half height
Gu1 = imgaussfilt3(u1,Gstd);
Gu2 = imgaussfilt3(u2,Gstd);
Gu3 = imgaussfilt3(u3,Gstd);

% SGS stress tensor = <uiuj> - <ui><uj>
% <.> = filtering operation
T11 = imgaussfilt3(u1.*u1,Gstd) - Gu1.*Gu1;
T12 = imgaussfilt3(u1.*u2,Gstd) - Gu1.*Gu2;
T13 = imgaussfilt3(u1.*u3,Gstd) - Gu1.*Gu3;
T22 = imgaussfilt3(u2.*u2,Gstd) - Gu2.*Gu2;
T23 = imgaussfilt3(u2.*u3,Gstd) - Gu2.*Gu3;
T33 = imgaussfilt3(u3.*u3,Gstd) - Gu3.*Gu3; % symmetric

% specify grid widths
[ny, nx, nz] = size(dat,'U');
hx = 8*pi/(nx-1);
hy = 2/(ny-1);
hz = 3*pi/(nz-1);

% calculate gradients of filtered data
[dUdy,dUdx,dUdz] = gradient(Gu1,hy,hx,hz);
[dVdy,dVdx,dVdz] = gradient(Gu2,hy,hx,hz);
[dWdy,dWdx,dWdz] = gradient(Gu3,hy,hx,hz);

% strain rate tensor
S11 = dUdx;
S12 = (dUdy + dVdx)/2;
S13 = (dUdz + dWdx)/2;
S22 = dVdy;
S23 = (dVdz + dWdy)/2;
S33 = dWdz; % symmetric

% calculate invariant
Snorm = sqrt(2*(S11.*S11 + 2*S12.*S12 + 2*S13.*S13 + S22.*S22 + 2*S23.*S23 + S33.*S33));

% ---- Smagorinsky Model ----
Delta = (hx*hy*hz)^(1/3);
cs = 0.16;

% stress tensor
repc = -2*(cs*Delta)^2;
Tsmag11 = repc*Snorm.*S11;
Tsmag12 = repc*Snorm.*S12;
Tsmag13 = repc*Snorm.*S13;
Tsmag22 = repc*Snorm.*S22;
Tsmag23 = repc*Snorm.*S23;
Tsmag33 = repc*Snorm.*S33;

% % if we want to load into variables
% x1 = dat.X;
% x2 = dat.Y;
% x3 = dat.Z;

% ---- Plotting ----
figure(1), clf
contourf(dat.X(:,:,1), dat.Y(:,:,1), T12(:,:,1), 30, 'LineStyle', 'none');
cbar = colorbar;
cbar.Limits = [-0.013 0.013];

figure(2), clf
contourf(dat.X(:,:,1), dat.Y(:,:,1), Tsmag12(:,:,1), 30, 'LineStyle', 'none');
cbar = colorbar;
cbar.Limits = [-0.013 0.013];