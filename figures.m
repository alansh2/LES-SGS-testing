% figures.m ...

clear all;
close all;

pc = parcluster('local');
parpool(pc, 3);

% load data
df1 = matfile('data.mat');
df2 = matfile('resolved.mat');
[ny,nx,nz] = size(df1,'U');
z0 = uint16(nz/2); % select center C/S
U = df1.U(:,:,z0);
Ures = df2.Gu1(:,:,z0);
Usgs = U - Ures;
info = {U, 'DNS'; Ures, 'resolved'; Usgs, 'unresolved'};

X = df1.X(:,:,z0);
Y = df1.Y(:,:,z0);

cblim = [min(U,[],'all') max(U,[],'all')];

parfor idx = 1:3
    [data,name] = info{idx,:};
    figure(idx), clf
    axes('Units', 'normalized', 'Position', [0 0 1 1])
    contourf(X, Y, data, 25000, 'LineStyle', 'none')
    colormap jet
    caxis(cblim)
    axis([X(1) X(1,nx) Y(1) Y(ny)])
    set(gca,'Visible','off')

    f = gcf;
    f.PaperUnits = 'points';
    f.PaperPosition = [0 0 921.6 172.8];
    print(name,'-dpng','-r300')
end
