function ev = Smagorinsky(Delta,S);
% Returns kinematic eddy viscosity
% calculate invariant
Snorm = sqrt(2*(S{1,1}.*S{1,1} + 2*S{1,2}.*S{1,2} + 2*S{1,3}.*S{1,3} ...
    + S{2,2}.*S{2,2} + 2*S{2,3}.*S{2,3} + S{3,3}.*S{3,3}));

cs = 0.16; % Smagorinsky coefficient
ev = -2*(cs*Delta)^2*Snorm; % eddy viscosity