% TensorPlots.m ... Plotting function for master.m
% fields: Tensor components requested to plot
%    X,Y: Coordinate grid
%  T,mod: Exact and model stress tensors
% zslice: Location of 2D plane taken for plotting
%  axlim: 2D Axis limits
function TensorPlots(fields,X,Y,T,mod,zslice,axlim,dependency,symmetry);
    name = mod.Name; % model name
    for idx = 1:numel(fields)
        % generate requested colorbar limits
        [lim1, lim2] = setColorbar(T,mod,fields{idx},zslice,dependency,symmetry);
        
        figure(idx); clf;
        
        subplot(2,1,1)
        contourf(X(:,:,zslice), Y(:,:,zslice), T.(fields{idx})(:,:,zslice), 30, 'LineStyle', 'none');
        set(gca, 'FontSize', 11);
        title([fields{idx} ' Exact'], 'FontSize', 13, 'FontWeight', 'bold');
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
        axis(axlim);
        set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');
        colorbar('FontSize', 12);
        caxis(lim1);
    
        subplot(2,1,2)
        contourf(X(:,:,zslice), Y(:,:,zslice), mod.(fields{idx})(:,:,zslice), 30, 'LineStyle', 'none');
        set(gca, 'FontSize', 11);
        title([fields{idx} ' ' name], 'FontSize', 13, 'FontWeight', 'bold');
        xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
        axis(axlim);
        set(gca, 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on');
        colorbar('FontSize', 12);
        caxis(lim2);

        fig = gcf;
        fig.PaperUnits = 'points';
        fig.PaperPosition = [0 0 1000 500];
        print([fields{idx} '_' name],'-dpng','-r300')
    end
end

function [lim1, lim2] = setColorbar(T,mod,Tij,zslice,dependency,symmetry);
    if strcmp(dependency,'match')
        up = max([max(T.(Tij)(:,:,zslice),[],'all'), max(mod.(Tij)(:,:,zslice),[],'all')]);
        lo = min([min(T.(Tij)(:,:,zslice),[],'all'), min(mod.(Tij)(:,:,zslice),[],'all')]);
        if strcmp(symmetry,'symmetric')
            fprintf('Using matching symmetric colorbar limits\n');
            bndy = max(abs([lo up]));
            lim1 = [-bndy bndy];
        elseif strcmp(symmetry,'none')
            fprintf('Using matching automatic colorbar limits\n');
            lim1 = [lo up];
        else
            fprintf('Unrecognized input argument "%s"\nUsing matching automatic colorbar limits\n',symmetry);
            lim1 = [lo up];
        end
        lim2 = lim1;
    elseif strcmp(dependency,'none')
        if strcmp(symmetry,'symmetric')
            fprintf('Using independent symmetric colorbar limits\n');
            bndy = max(abs([min(T.(Tij)(:,:,zslice),[],'all'), max(T.(Tij)(:,:,zslice),[],'all')]));
            lim1 = [-bndy bndy];
            bndy = max(abs([min(mod.(Tij)(:,:,zslice),[],'all'), max(mod.(Tij)(:,:,zslice),[],'all')]));
            lim2 = [-bndy bndy];
        elseif strcmp(symmetry,'none')
            fprintf('Using automatic colorbar limits\n');
            lim1 = 'auto'; lim2 = 'auto';
        else
            fprintf('Unrecognized input argument "%s"\nUsing automatic colorbar limits\n',dependency);
            lim1 = 'auto'; lim2 = 'auto';
        end
    end
end