% ExactSGS.m ... Part of a priori test
% Compute exact SGS stress tensor and strain rate tensor
function [T,S] = ExactSGS(u1,u2,u3,std,hx,hy,hz);
% Gaussian filter
Gu1 = imgaussfilt3(u1,std);
Gu2 = imgaussfilt3(u2,std);
Gu3 = imgaussfilt3(u3,std);
% SGS stress:      ___    _ _
%            Tij = uiuj = uiuj
T = cell(3,3);
T{1,1} = imgaussfilt3(u1.*u1,std) - Gu1.*Gu1;
T{1,2} = imgaussfilt3(u1.*u2,std) - Gu1.*Gu2;
T{1,3} = imgaussfilt3(u1.*u3,std) - Gu1.*Gu3;
T{2,2} = imgaussfilt3(u2.*u2,std) - Gu2.*Gu2;
T{2,3} = imgaussfilt3(u2.*u3,std) - Gu2.*Gu3;
T{3,3} = imgaussfilt3(u3.*u3,std) - Gu3.*Gu3;

% get gradients for S
[dUdy,dUdx,dUdz] = gradient(Gu1,hy,hx,hz);
clear Gu1;
[dVdy,dVdx,dVdz] = gradient(Gu2,hy,hx,hz);
clear Gu2;
[dWdy,dWdx,dWdz] = gradient(Gu3,hy,hx,hz);
clear Gu3;

% strain rate tensor
S = cell(3,3);
S{1,1} = dUdx; clear dUdx;
S{1,2} = 0.5*(dUdy + dVdx); clear dUdy dVdx;
S{1,3} = 0.5*(dUdz + dWdx); clear dUdz dWdx;
S{2,2} = dVdy; clear dVdy;
S{2,3} = 0.5*(dVdz + dWdy); clear dVdz dWdy;
S{3,3} = dWdz; clear dWdz;