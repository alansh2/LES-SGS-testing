% SGSmodels.m ... 

clear all;
close all;

load('properties.mat'); % contains integral length scale
df1 = matfile('data.mat');
df2 = importdata('profiles.txt',' ',2);

X = df1.X;
Y = df1.Y;
x = X(1,:,1);
y = Y(:,1,1);
ypl = df2.data(1:end-1,1);

T = matfile('T.mat');
S = matfile('S.mat');
g = matfile('gradvel.mat');

Smag = @() Smagorinsky(Delta,S,y,ypl);
WALE = @() WallAdapting(Delta,S,g);

mod = closure(Smag,S);

%% Plots
% setup
nz = size(df1,'X',3);
zslice = uint16(nz/2);
axlim = [min(x) max(x) min(y) max(y)];
fields = {'T12','T13','T23'}; % choose components
dependency = 'match'; % match, none
symmetry = 'symmetric'; % symmetric, none
% get plots
parpool(parcluster('local'), numel(fields));
TensorPlots(fields,X,Y,T,mod,zslice,axlim,dependency,symmetry);