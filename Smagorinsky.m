function ev = Smagorinsky(Delta,S);
    % Returns kinematic eddy viscosity
    % calculate invariant
    Snorm = sqrt(2*(S.S11.*S.S11 + 2*S.S12.*S.S12 + 2*S.S13.*S.S13 ...
        + S.S22.*S.S22 + 2*S.S23.*S.S23 + S.S33.*S.S33));

    cs = 0.16; % Smagorinsky coefficient
    ev = -2*(cs*Delta)^2*Snorm; % eddy viscosity
end