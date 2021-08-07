% ChannelResolvent.m ... 

clear
close all

% read data
df = matfile('data.mat');
[ny,nx,nz] = size(df,'X');
y = df.Y(:,1,1)';
dx = df.X(1,2,1) - df.X(1,1,1);
dz = df.Z(1,1,2) - df.Z(1,1,1);
h = y(end);
L = df.X(1,nx,1);
D = df.Z(1,1,nz);

% SGS parameters
P = 4*h + 2*D; % wetted perimeter
d_h = 8*h*D/P; % hydraulic diameter
% characteristic length scale as defined by
% https://www.cfd-online.com/Wiki/Hydraulic_diameter
Delta = 0.07*d_h;
Gstd = sqrt(Delta^2/12)*[1/dx 1 1/dz]; % standard deviation of Gaussian filter
save('properties.mat','Delta','Gstd');

% resolved scales from Gaussian filter
Gu1 = GaussianFilter(y,df.U,Gstd);
Gu2 = GaussianFilter(y,df.V,Gstd);
Gu3 = GaussianFilter(y,df.W,Gstd);
save('resolved.mat','Gu1','Gu2','Gu3','-v7.3');
