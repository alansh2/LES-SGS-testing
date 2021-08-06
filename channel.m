% channel.m ... Get and save channel DNS data from JHTDB
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
% Turbmat tools from John Hopkins Turbulence Databases
%

clear all;
close all;

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
pc = parcluster('local');
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));

%% Properties
% snapshot time
time = 0.0;

% channel geometry
h = 1; % channel half height
L = 8*pi; % streamwise length
D = 3*pi; % channel depth

% construct grid
% build y from imported y+
df = importdata('profiles.txt',' ',2);
ypl = df.data(1:256,1)';
v = 5e-5; ut = 4.9968e-2;
yL = ypl*v/ut - h;
y = [yL -flip(yL)];

nx = 2048;
ny = length(y);
nz = 1536;
npoints = nx*ny*nz;
nmax = 4096;
n = npoints/nmax;

x = linspace(0, L, nx);
z = linspace(0, D, nz);
[X, Y, Z] = meshgrid(x, y, z);

% adjust point structure
points = zeros(3, nmax, n);
points(1,:,:) = reshape(X(:), nmax, n);
points(2,:,:) = reshape(Y(:), nmax, n);
points(3,:,:) = reshape(Z(:), nmax, n);

% point sampling interpolation parameters
SInt = 'Lag4';
TInt = 'None';

%% Get Velocities
% parGetVel.m ...
fprintf('\nRequesting velocity at %is...\n',time);
authkey = 'edu.jhu.pha.turbulence.testing-201406';
% initialize velocities
u = zeros(nmax,n);
v = zeros(nmax,n);
w = zeros(nmax,n);
% get velocity at each point
parfor i = 1:n
    pslice = points(:,:,i);
    vel = getVelocity(authkey, 'channel', time, SInt, TInt, nmax, pslice);
    u(:,i) = vel(1,:);
    v(:,i) = vel(2,:);
    w(:,i) = vel(3,:);
end
clear points

U = reshape(u(:), ny, nx, nz);
clear u;
V = reshape(v(:), ny, nx, nz);
clear v;
W = reshape(w(:), ny, nx, nz);
clear w;

save('data.mat','U','V','W','X','Y','Z','-v7.3');
