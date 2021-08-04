%% Master program for data collection and analysis
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
mkdir temp

clear all;
close all;

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
pc = parcluster('local');
parpool(pc, 4);%str2num(getenv('SLURM_CPUS_ON_NODE')));

%% Properties
time = 0.0;

load('y.mat');
nx = 128;%2048;
ny = length(y); ny = 32;
nz = 96;%1536;
npoints = nx*ny*nz;
nmax = 4096;
n = npoints/nmax;

% construct grid
h = 1; % channel half height
L = 8*pi; % streamwise length
D = 3*pi; % channel depth
x = linspace(0, L, nx);
% y-grid from y.mat
y = linspace(-h, h, ny);
z = linspace(0, D, nz);
[X, Y, Z] = meshgrid(x, y, z);

dx = x(2) - x(1);
dz = z(2) - z(1);

% adjust point structure
points = zeros(3, nmax, n);
points(1,:,:) = reshape(X(:), nmax, n);
points(2,:,:) = reshape(Y(:), nmax, n);
points(3,:,:) = reshape(Z(:), nmax, n);

% point sampling interpolation parameters
SInt = 'Lag4';
TInt = 'None';

% SGS parameters
P = 4*h + 2*D; % wetted perimeter
dh = 8*h*D/P; % hydraulic diameter
% characteristic length scale 'Delta'
% as defined by https://www.cfd-online.com/Wiki/Hydraulic_diameter
Delta = 0.07*dh;
Gstd = sqrt(Delta^2/12)*[1/dx 1 1/dz]; % standard deviation of Gaussian filter

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

%% Exact Solution
[T,S,g] = ExactSGS(U,V,W,Gstd,x,y,z);
clear U V W;

correlation = @(tau) tau.T11.*S.S11 + tau.T22.*S.S22 + tau.T33.*S.S33 ...
    + 2*(tau.T12.*S.S12 + tau.T13.*S.S13 + tau.T23.*S.S23);

%% SGS Models
Smag = @() Smagorinsky(Delta,S);
WALE = @() WallAdapting(Delta,S,g);

% eddy-viscosity closure
mod = closure(Smag,S);

% calculate angle between tensors
AB_IP = T.T11.*mod.T11 + 2*T.T12.*mod.T12 + 2*T.T13.*mod.T13 ...
    + T.T22.*mod.T22 + 2*T.T23.*mod.T23 + T.T33.*mod.T33;
Amag = sqrt(T.T11.^2 + 2*T.T12.^2 + 2*T.T13.^2 + T.T22.^2 + 2*T.T23.^2 + T.T33.^2);
Bmag = sqrt(mod.T11.^2 + 2*mod.T12.^2 + 2*mod.T13.^2 + mod.T22.^2 + 2*mod.T23.^2 + mod.T33.^2);
angle = acosd(AB_IP./Amag./Bmag);
disp(mean(angle,'all'));

%% Plots
% setup
zslice = uint16(nz/2);
axlim = [min(x) max(x) min(y) max(y)];
fields = {'T12','T13','T23'}; % choose components
dependency = 'match'; % match, none
symmetry = 'symmetric'; % symmetric, none
% get plots
TensorPlots(fields,X,Y,T,mod,zslice,axlim,dependency,symmetry);
