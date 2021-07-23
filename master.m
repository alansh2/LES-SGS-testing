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
mkdir temp

clear all;
close all;

% start the matlabpool with maximum available workers
% control how many workers by setting ntasks in your sbatch script
pc = parcluster('local');
parpool(pc, 4);%str2num(getenv('SLURM_CPUS_ON_NODE')));

%% Properties
t0 = 12.870; % starting time
dt = 0.0013;
Nt = 1; % number of snapshots to average

nx = 128;%2048;
ny = 32;%512;
nz = 96;%1536;
npoints = nx*ny*nz;
nmax = 4096;
n = npoints/nmax;

% construct grid
h = 1; % channel half height
L = 8*pi; % streamwise length
D = 3*pi; % channel depth
x = linspace(0, L, nx);
y = linspace(-h, h, ny);
z = linspace(0, D, nz);
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

% point sampling interpolation parameters
Spatial = 'Lag4';
Temporal = 'None';

% SGS parameters
P = 4*h + 2*D; % wetted perimeter
Delta = 2*h*D/P; % characteristic length scale
Gstd = sqrt(Delta^2/12)*[1/dx 1/dy 1/dz]; % standard deviation of Gaussian filter

%% Core Task
% purges exhausted variables
exact = zeros(ny, nx, nz);
model = zeros(ny, nx, nz);
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
    [T,S,g] = ExactSGS(U,V,W,Gstd,dx,dy,dz);
    clear U V W;
    
    corr = T.T11.*S.S11 + 2*T.T12.*S.S12 + 2*T.T13.*S.S13 ...
        + T.T22.*S.S22 + 2*T.T23.*S.S23 + T.T33.*S.S33;
    
    exact = exact + corr;
    
    % ---- SGS Models ----
    Smag = @() Smagorinsky(Delta,S);
    % WALE = @() WallAdaptingLE(Delta,S,g); % under construction
    
    % eddy-viscosity closure
    mod = closure(Smag,S);
    corr = mod.T11.*S.S11 + 2*mod.T12.*S.S12 + 2*mod.T13.*S.S13 ...
        + mod.T22.*S.S22 + 2*mod.T23.*S.S23 + mod.T33.*S.S33;
    
    model = model + corr;
    
    
%     % calculate angle between tensors
%     AB_IP = T.T11.*mod.T11 + 2*T.T12.*mod.T12 + 2*T.T13.*mod.T13 ...
%         + T.T22.*mod.T22 + 2*T.T23.*mod.T23 + T.T33.*mod.T33;
%     Amag = sqrt(T.T11.^2 + 2*T.T12.^2 + 2*T.T13.^2 + T.T22.^2 + 2*T.T23.^2 + T.T33.^2);
%     Bmag = sqrt(mod.T11.^2 + 2*mod.T12.^2 + 2*mod.T13.^2 + mod.T22.^2 + 2*mod.T23.^2 + mod.T33.^2);
%     angle = acosd(AB_IP./Amag./Bmag);
%     disp(mean(angle,'all'));
end

%% Plots
% setup colorbar limits
zi = uint16(nz/2);
up = max([max(T.T13(:,:,zi),[],'all'), max(mod.T13(:,:,zi),[],'all')]);
lo = min([min(T.T13(:,:,zi),[],'all'), min(mod.T13(:,:,zi),[],'all')]);

% can plot other components of stress tensor
% important that the components are of the same order of magnitude
% use caxis to compare on same levels
subplot(2,1,1)
contourf(X(:,:,zi), Y(:,:,zi), T.T13(:,:,zi), 30, 'LineStyle', 'none');
set(gca, 'FontSize', 11);
title('T13 Exact', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
colorbar('FontSize', 12);
% caxis([lo up]);
axis([min(x) max(x) min(y) max(y)]);
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');

subplot(2,1,2)
contourf(X(:,:,zi), Y(:,:,zi), mod.T13(:,:,zi), 30, 'LineStyle', 'none');
set(gca, 'FontSize', 11);
title('T13 Model', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
colorbar('FontSize', 12);
% caxis([lo up]);
axis([min(x) max(x) min(y) max(y)]);
set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');

fig = gcf;
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 1000 500];
print('T13','-dpng','-r0')
