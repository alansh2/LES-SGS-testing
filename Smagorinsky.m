function ev = Smagorinsky(Delta,S);
% Returns kinematic eddy viscosity
% calculate invariant
Snorm = sqrt(2*(S.S11.^2 + S.S22.^2 + S.S33.^2 + 2*(S.S12.^2 + S.S13.^2 + S.S23.^2)));

cs = 0.18; % Smagorinsky coefficient
ev = (cs*Delta)^2*Snorm; % eddy viscosity
