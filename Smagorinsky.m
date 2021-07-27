function ev = Smagorinsky(Delta,S);
    % Returns kinematic eddy viscosity
    % calculate invariant
    Snorm = sqrt(2*(S.S11.^2 + 2*S.S12.^2 + 2*S.S13.^2 ...
        + S.S22.^2 + 2*S.S23.^2 + S.S33.^2));

    cs = 0.16; % Smagorinsky coefficient
    ev = -2*(cs*Delta)^2*Snorm; % eddy viscosity
end