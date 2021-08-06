function ev = Smagorinsky(Delta,S,y,ypl);

% An eddy-viscosity closure model. Returns eddy viscosity.

% calculate invariant
Snorm = sqrt(2*(S.S11.^2 + S.S22.^2 + S.S33.^2 + 2*(S.S12.^2 + S.S13.^2 + S.S23.^2)));

cs = 0.17; % Smagorinsky coefficient
kappa = 0.41;
hy = -abs(y) + y(end);
VDS = -exp(-[ypl; flip(ypl)]/26) + 1;
l = min(kappa*hy.*VDS,cs*Delta);
ev = l.^2.*Snorm; % eddy viscosity
