% strain.m ... 

clear
close all

load('resolved.mat')
df = matfile('data.mat');
[ny,nx,nz] = size(df,'X');
x = df.X(1,:,1);
y = df.Y(:,1,1)';
z = reshape(df.Z(1,1,:),1,nz);

% velocity gradients
g = matfile('gradvel.mat','Writable',true);
[g.dUdx,g.dUdy,g.dUdz] = gradient(Gu1,x,y,z);
clear Gu1;
[g.dVdx,g.dVdy,g.dVdz] = gradient(Gu2,x,y,z);
clear Gu2;
[g.dWdx,g.dWdy,g.dWdz] = gradient(Gu3,x,y,z);
clear Gu3;

% strain rate tensor
Sbar = matfile('S.mat','Writable',true);
Sbar.S11 = g.dUdx;
Sbar.S12 = 0.5*(g.dUdy + g.dVdx);
Sbar.S13 = 0.5*(g.dUdz + g.dWdx);
Sbar.S22 = g.dVdy;
Sbar.S23 = 0.5*(g.dVdz + g.dWdy);
Sbar.S33 = g.dWdz;
