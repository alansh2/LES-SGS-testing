% ExactSGS.m ... Part of a priori test
% Compute exact SGS stress tensor and strain rate tensor
function [tau,Sbar,g] = ExactSGS(u1,u2,u3,std,x,y,z);
% Gaussian filter
Gu1 = GaussianFilter(y,u1,std);
Gu2 = GaussianFilter(y,u2,std);
Gu3 = GaussianFilter(y,u3,std);
% SGS stress:   ___    _ _
%         Tij = uiuj = uiuj
tau = matfile('./temp/T.mat','Writable',true);
tau.T11 = GaussianFilter(y,u1.*u1,std) - Gu1.*Gu1;
tau.T12 = GaussianFilter(y,u1.*u2,std) - Gu1.*Gu2;
tau.T13 = GaussianFilter(y,u1.*u3,std) - Gu1.*Gu3;
tau.T22 = GaussianFilter(y,u2.*u2,std) - Gu2.*Gu2;
tau.T23 = GaussianFilter(y,u2.*u3,std) - Gu2.*Gu3;
tau.T33 = GaussianFilter(y,u3.*u3,std) - Gu3.*Gu3;

% get gradients for 'Sbar'
% velocity gradient tensor 'g'
g = matfile('./temp/gradvel.mat','Writable',true);
[g.dUdx,g.dUdy,g.dUdz] = gradient(Gu1,x,y,z);
clear Gu1;
[g.dVdx,g.dVdy,g.dVdz] = gradient(Gu2,x,y,z);
clear Gu2;
[g.dWdx,g.dWdy,g.dWdz] = gradient(Gu3,x,y,z);
clear Gu3;

% strain rate tensor
Sbar = matfile('./temp/S.mat','Writable',true);
Sbar.S11 = g.dUdx;
Sbar.S12 = 0.5*(g.dUdy + g.dVdx);
Sbar.S13 = 0.5*(g.dUdz + g.dWdx);
Sbar.S22 = g.dVdy;
Sbar.S23 = 0.5*(g.dVdz + g.dWdy);
Sbar.S33 = g.dWdz;