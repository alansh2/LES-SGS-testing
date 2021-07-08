% testLES_SGS.m - LES SGS Model Test From DNS
% 
% 

close all
clear all

% load DNS
dat = matfile('data.mat');

% SGS stress tensor = <uiuj> - <ui><uj>
% <.> = filtering operation
% Tij = 

% Gaussian filter
Gstd = 0.5; % half channel half height
Gu1 = imgaussfilt3(dat.U,Gstd);
Gu2 = imgaussfilt3(dat.V,Gstd);
Gu3 = imgaussfilt3(dat.W,Gstd);

% specify grid widths
[ny, nx, nz] = size(dat,'U');
hx = 8*pi/(nx-1);
hy = 2/(ny-1);
hz = 3*pi/(nz-1);
% size(dat.VAR) = ny x nx x nz
[dUdy,dUdx,dUdz] = gradient(Gu1,hy,hx,hz);
[dVdy,dVdx,dVdz] = gradient(Gu2,hy,hx,hz);
[dWdy,dWdx,dWdz] = gradient(Gu3,hy,hx,hz);

S11 = dUdx;
S12 = (dUdy + dVdx)/2;
S13 = (dUdz + dWdx)/2;
S22 = dVdy;
S23 = (dVdz + dWdy)/2;
S33 = dWdz;
% symmetric

Snorm = sqrt(2*(S11.*S11 + 2*S12.*S12 + 2*S13.*S13 + S22.*S22 + 2*S23.*S23 + S33.*S33));

Delta = (hx*hy*hz)^(1/3);
cs = 0.16;

% stress tensor
repc = -2*(cs*Delta)^2;
T11 = repc*Snorm.*S11;
T12 = repc*Snorm.*S12;
T13 = repc*Snorm.*S13;
T22 = repc*Snorm.*S22;
T23 = repc*Snorm.*S23;
T33 = repc*Snorm.*S33;

% plotting
figure(1)
contourf(dat.X(:,:,1), dat.Y(:,:,1), T12(:,:,1), 30, 'LineStyle', 'none');
colorbar
