% master.m ... Master program for data collection and analysis
%
% Written by:
%
% Alan Hong
% University of Illinois - Urbana Champaign
% Saxton-Fox Group @ UIUC
%

mkdir data

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
x = linspace(0, 8*pi, nx);
y = linspace(-1, 1, ny);
z = linspace(0, 3*pi, nz);
[X, Y, Z] = meshgrid(x, y, z);

% specify grid spacing
hx = x(2) - x(1);
hy = y(2) - y(1);
hz = z(2) - z(1);

% adjust point structure
points = zeros(3, nmax, n);
points(1,:,:) = reshape(X(:), nmax, n);
points(2,:,:) = reshape(Y(:), nmax, n);
points(3,:,:) = reshape(Z(:), nmax, n);

% interpolation parameters
Spatial = 'Lag4';
Temporal = 'None';

% get velocities at each point for Nt snapshots in time
for idx = 0:Nt-1
    [u, v, w] = parGetVel(t0+5*dt*idx, nmax, n, points, Spatial, Temporal);
    U = reshape(u(:), ny, nx, nz);
    clear u;
    V = reshape(v(:), ny, nx, nz);
    clear v;
    W = reshape(w(:), ny, nx, nz);
    clear w;
    save(['./data/t_' num2str(idx) '.mat'],'U','V','W','-v7.3');
    tau = Smagorinsky(U, V, W, hx, hy, hz);
end
