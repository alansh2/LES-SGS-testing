% ExactSGS.m ... Part of a priori test
% Compute exact SGS stress tensor and strain rate tensor
function [tau,Sbar,g] = ExactSGS(u1,u2,u3,std,hx,hy,hz);
    % Gaussian filter
    Gu1 = imgaussfilt3(u1,std);
    Gu2 = imgaussfilt3(u2,std);
    Gu3 = imgaussfilt3(u3,std);
    % SGS stress:   ___    _ _
    %         Tij = uiuj = uiuj
    tau = matfile('./temp/T.mat','Writable',true);
    tau.T11 = imgaussfilt3(u1.*u1,std) - Gu1.*Gu1;
    tau.T12 = imgaussfilt3(u1.*u2,std) - Gu1.*Gu2;
    tau.T13 = imgaussfilt3(u1.*u3,std) - Gu1.*Gu3;
    tau.T22 = imgaussfilt3(u2.*u2,std) - Gu2.*Gu2;
    tau.T23 = imgaussfilt3(u2.*u3,std) - Gu2.*Gu3;
    tau.T33 = imgaussfilt3(u3.*u3,std) - Gu3.*Gu3;
    
    % get gradients for 'Sbar'
    % velocity gradient tensor 'g'
    g = matfile('./temp/gradvel.mat','Writable',true);
    [g.dUdy,g.dUdx,g.dUdz] = gradient(Gu1,hy,hx,hz);
    clear Gu1;
    [g.dVdy,g.dVdx,g.dVdz] = gradient(Gu2,hy,hx,hz);
    clear Gu2;
    [g.dWdy,g.dWdx,g.dWdz] = gradient(Gu3,hy,hx,hz);
    clear Gu3;

    % strain rate tensor
    Sbar = matfile('./temp/S.mat','Writable',true);
    Sbar.S11 = g.dUdx;
    Sbar.S12 = 0.5*(g.dUdy + g.dVdx);
    Sbar.S13 = 0.5*(g.dUdz + g.dWdx);
    Sbar.S22 = g.dVdy;
    Sbar.S23 = 0.5*(g.dVdz + g.dWdy);
    Sbar.S33 = g.dWdz;
end