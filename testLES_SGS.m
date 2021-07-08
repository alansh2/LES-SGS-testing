% testLES_SGS.m - LES SGS Model Test From DNS
%
%     y
%     |  ,----------------------,
%     | /|                     /|
%     |/ |                    / |
%     2--+-------------------:  |
%     |  | z                 |  |
%     |  |/                  |  |
% 512 | 3*pi-----------------+--;
%     | /                    | / 1536
%     |/                     |/
%     +---------------------8*pi--->x
%               2048
% 

close all
clear all

% load DNS
dat = matfile('data.mat');

% suppose we have 'u' extracted [U; V; W]
% assume single frame data is temporal average
% SGS stress tensor = <uiuj> - <ui><uj>
% <.> = filtering operation

Gstd = 0.5; % half channel half height
u1 = imgaussfilt3(dat.U,Gstd);
u2 = imgaussfilt3(dat.V,Gstd);
u3 = imgaussfilt3(dat.W,Gstd);

% specify grid widths
[ny, nx, nz] = size(dat,'U');
hx = 8*pi/(nx-1);
hy = 2/(ny-1);
hz = 3*pi/(nz-1);
% size(dat.VAR) = ny x nx x nz
[dUdy,dUdx,dUdz] = gradient(u1,hy,hx,hz);
[dVdy,dVdx,dVdz] = gradient(u2,hy,hx,hz);
[dWdy,dWdx,dWdz] = gradient(u3,hy,hx,hz);

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
T11 = -2*(cs*Delta)^2*Snorm.*S11;
T12 = -2*(cs*Delta)^2*Snorm.*S12;
T13 = -2*(cs*Delta)^2*Snorm.*S13;
T22 = -2*(cs*Delta)^2*Snorm.*S22;
T23 = -2*(cs*Delta)^2*Snorm.*S23;
T33 = -2*(cs*Delta)^2*Snorm.*S33;
