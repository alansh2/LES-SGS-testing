% parGetVel.m ... Parallel Velocity Sampling Algorithm
% 
% Accelerated data sampling from JHU Turbulence Database Cluster
%
% Algorithm by:
%
% Alan Hong
% University of Illinois at Urbana-Champaign
% Saxton-Fox Group @ UIUC
%
% Ricky Hsu
% Arizona State University
%
% Turbmat tools from John Hopkins University Turbulence Databases
%
function [u, v, w] = parGetVel(time, nmax, n, points, SInt, TInt);
fprintf('\nRequesting velocity at T+%is...\n',time);
authkey = 'edu.jhu.pha.turbulence.testing-201406';
% initialize velocities
u = zeros(nmax,n);
v = zeros(nmax,n);
w = zeros(nmax,n);
parfor i = 1:n
    pslice = points(:,:,i);
    vel = getVelocity(authkey, 'channel', time, SInt, TInt, nmax, pslice);
    u(:,i) = vel(1,:);
    v(:,i) = vel(2,:);
    w(:,i) = vel(3,:);
end