% Parallel Compute Algorithm
%
% Accelerated data sampling from JHU Turbulence Database Cluster
%
% Algorithm by:
%
% Alan Hong
% University of Illinois - Urbana Champaign
% Saxton-Fox Group @ UIUC
%
% Ricky Hsu
% Arizona State University
%
% Turbmat from John Hopkins University Turbulence Databases
%

clear all;
close all;

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
pc = parcluster('local');
parpool(pc, 4); % str2num(getenv('SLURM_CPUS_ON_NODE')));

% Turbmat tools
authkey = 'edu.jhu.pha.turbulence.testing-201406';
dataset = 'channel';

% ---- Temporal Interpolation Options ----
NoTInt   = 'None' ; % No temporal interpolation
PCHIPInt = 'PCHIP'; % Piecewise cubic Hermit interpolation in time

% ---- Spatial Interpolation Flags for getVelocity & getVelocityAndPressure ----
NoSInt = 'None'; % No spatial interpolation
Lag4   = 'Lag4'; % 4th order Lagrangian interpolation in space
Lag6   = 'Lag6'; % 6th order Lagrangian interpolation in space
Lag8   = 'Lag8'; % 8th order Lagrangian interpolation in space

% ---- Spatial Differentiation & Interpolation Flags for getVelocityGradient & getPressureGradient ----
FD4NoInt = 'None_Fd4' ; % 4th order finite differential scheme for grid values, no spatial interpolation
FD6NoInt = 'None_Fd6' ; % 6th order finite differential scheme for grid values, no spatial interpolation
FD8NoInt = 'None_Fd8' ; % 8th order finite differential scheme for grid values, no spatial interpolation
FD4Lag4  = 'Fd4Lag4'  ; % 4th order finite differential scheme for grid values, 4th order Lagrangian interpolation in space

% ---- Spline interpolation and differentiation Flags for getVelocity,
% getPressure, getVelocityGradient, getPressureGradient,
% getVelocityHessian, getPressureHessian
M1Q4   = 'M1Q4'; % Splines with smoothness 1 (3rd order) over 4 data points. Not applicable for Hessian.
M2Q8   = 'M2Q8'; % Splines with smoothness 2 (5th order) over 8 data points.
M2Q14   = 'M2Q14'; % Splines with smoothness 2 (5th order) over 14 data points.

time = 0;
dt = 0.0013;
Nt = 100;

nx = 128; % 2048;
ny = 32; % 512;
nz = 8; % 1536;
npoints = nx*ny*nz;
nmax = 4096;
n = npoints/nmax;
points = zeros(3, nmax, n);

x = linspace(0, 8*pi, nx);
y = linspace(-1, 1, ny);
z = linspace(0, 3*pi, nz);
[X, Y, Z] = meshgrid(x, y, z);
points(1,:,:) = reshape(X(:), nmax, n);
points(2,:,:) = reshape(Y(:), nmax, n);
points(3,:,:) = reshape(Z(:), nmax, n);

% initialize velocities
u = zeros(nmax,n);
v = zeros(nmax,n);
w = zeros(nmax,n);

% get velocities at each point
fprintf('\nRequesting velocity at (%ix%ix%i) points...\n',nx,ny,nz);
parfor i = 1:n
    pslice = points(:,:,i);
    for j = 0:Nt-1
        vel = getVelocity(authkey, dataset, time+5*dt*j, Lag4, NoTInt, nmax, pslice);
        u(:,i) = u(:,i) + vel(1,:)';
        v(:,i) = v(:,i) + vel(2,:)';
        w(:,i) = w(:,i) + vel(3,:)';
    end
end

clear points;

% package velocities
U = reshape(u(:), ny, nx, nz)/Nt;
V = reshape(v(:), ny, nx, nz)/Nt;
W = reshape(w(:), ny, nx, nz)/Nt;

save('data.mat','U','V','W','X','Y','Z','-v7.3')
