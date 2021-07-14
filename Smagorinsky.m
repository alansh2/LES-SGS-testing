function tau = Smagorinsky(u1, u2, u3, hx, hy, hz);
% pass through Gaussian filter then calculate gradient
% clear obsolete variables to limit memory usage
Gstd = 0.5;
Gu1 = imgaussfilt3(u1,Gstd);
clear u1;
[dUdy,dUdx,dUdz] = gradient(Gu1,hy,hx,hz);
clear Gu1;
Gu2 = imgaussfilt3(u2,Gstd);
clear u2;
[dVdy,dVdx,dVdz] = gradient(Gu2,hy,hx,hz);
clear Gu2;
Gu3 = imgaussfilt3(u3,Gstd);
clear u3;
[dWdy,dWdx,dWdz] = gradient(Gu3,hy,hx,hz);
clear Gu3;

% strain rate tensor
S11 = dUdx;
clear dUdx;
S12 = 0.5*(dUdy + dVdx);
clear dUdy dVdx;
S13 = 0.5*(dUdz + dWdx);
clear dUdz dWdx;
S22 = dVdy;
clear dVdy;
S23 = 0.5*(dVdz + dWdy);
clear dVdz dWdy;
S33 = dWdz;
clear dWdz;

% calculate invariant
Snorm = sqrt(2*(S11.*S11 + 2*S12.*S12 + 2*S13.*S13 + S22.*S22 + 2*S23.*S23 + S33.*S33));

Delta = (hx*hy*hz)^(1/3);
cs = 0.16;

% stress tensor
nu = -2*(cs*Delta)^2*Snorm; % eddy viscosity
clear Snorm
tau = struct('T11',[],'T12',[],'T13',[],'T22',[],'T23',[],'T33',[]);
tau.T11 = nu.*S11; % clear S11;
tau.T12 = nu.*S12; % clear S12;
tau.T13 = nu.*S13; % clear S13;
tau.T22 = nu.*S22; % clear S22;
tau.T23 = nu.*S23; % clear S23;
tau.T33 = nu.*S33; % clear S33;