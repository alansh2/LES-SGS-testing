% master.m ... Master program for data collection and analysis
%
% Written by:
%
% Alan Hong
% University of Illinois at Urbana-Champaign
% Saxton-Fox Group @ UIUC
%
%
%
%
%
%
%
%

clear all;
close all;

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
pc = parcluster('local');
parpool(pc, 4); % str2num(getenv('SLURM_CPUS_ON_NODE')));

t0 = 0.0;
dt = 0.0013;
Nt = 2;

% properties
nx = 64; % 2048;
ny = 16; % 512;
nz = 48; % 1536;
npoints = nx*ny*nz;
nmax = 4096;
n = npoints/nmax;

% construct grid
h = 1; % channel half height
L = 8*pi; % streamwise length
W = 3*pi; % channel width
x = linspace(0, L, nx);
y = linspace(-h, h, ny);
z = linspace(0, W, nz);
[X, Y, Z] = meshgrid(x, y, z);

% specify grid spacing
dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);

% adjust point structure
points = zeros(3, nmax, n);
points(1,:,:) = reshape(X(:), nmax, n);
points(2,:,:) = reshape(Y(:), nmax, n);
points(3,:,:) = reshape(Z(:), nmax, n);

% interpolation parameters
Spatial = 'Lag4';
Temporal = 'None';

L_T = 0.07*W; % turbulent length scale
Gstd = L_T/6; % standard deviation of Gaussian filter

% ---- Core Task ----
% purge exhausted variables
for idx = 0:Nt-1
    % get velocity at each point
    [u, v, w] = parGetVel(t0+5*dt*idx, nmax, n, points, Spatial, Temporal);
    U = reshape(u(:), ny, nx, nz);
    clear u;
    V = reshape(v(:), ny, nx, nz);
    clear v;
    W = reshape(w(:), ny, nx, nz);
    clear w;
    
    % ---- A Priori Analysis ----
    [T,Sbar] = ExactSGS(U,V,W,Gstd,dx,dy,dz);
    
    % ---- Smagorinsky ----
    nu = Smagorinsky(L_T,Sbar);
    
    % eddy-viscosity closure
    T11 = nu.*Sbar{1,1}; Sbar(1,1) = {[]};
    T12 = nu.*Sbar{1,2}; Sbar(1,2) = {[]};
    T13 = nu.*Sbar{1,3}; Sbar(1,3) = {[]};
    T22 = nu.*Sbar{2,2}; Sbar(2,2) = {[]};
    T23 = nu.*Sbar{2,3}; Sbar(2,3) = {[]};
    T33 = nu.*Sbar{3,3}; Sbar(3,3) = {[]};
end