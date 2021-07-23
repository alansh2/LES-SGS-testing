% closure.m ... Eddy-Viscosity Closure
% Returns stress tensor from SGS model
function m = closure(model,S);
    nu = model(); % get eddy viscosity
    m = matfile('./temp/T_mod.mat','Writable',true);
    m.T11 = nu.*S.S11;
    m.T12 = nu.*S.S12;
    m.T13 = nu.*S.S13;
    m.T22 = nu.*S.S22;
    m.T23 = nu.*S.S23;
    m.T33 = nu.*S.S33;
end