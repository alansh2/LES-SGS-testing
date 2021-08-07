% stress.m ... 

clear
close all

% read data
load('properties.mat');
load('resolved.mat');
df = matfile('data.mat');
U = df.U;
V = df.V;
W = df.W;
y = df.Y(:,1,1)';

% SGS stress:   ___    _ _
%         Tij = uiuj = uiuj
tau = matfile('T.mat','Writable',true);
tau.T11 = GaussianFilter(y,U.*U,Gstd) - Gu1.*Gu1;
tau.T12 = GaussianFilter(y,U.*V,Gstd) - Gu1.*Gu2;
tau.T13 = GaussianFilter(y,U.*W,Gstd) - Gu1.*Gu3;
tau.T22 = GaussianFilter(y,V.*V,Gstd) - Gu2.*Gu2;
tau.T23 = GaussianFilter(y,V.*W,Gstd) - Gu2.*Gu3;
tau.T33 = GaussianFilter(y,W.*W,Gstd) - Gu3.*Gu3;
