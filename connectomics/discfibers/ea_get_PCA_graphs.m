function ea_get_PCA_graphs(data, var_names)
    % This function runs the PCA and plots diagnostic graphs that help to
    % interpret the PCs and to choose how many of them to retain
    % Data : a matrix with the data used, one variable per column, one line per
    % observation
    % var_names : names of the variables - cell array, 1 x N variables 
    
    % Components to consider for further plotting - avoids having too much plots
    if size(data, 2) <= 5 
        n2plot = size(data, 2);
    else
        n2plot = max([5, round(size(data, 2)./3)]);
    end 
    
    % Standardize data - z-score
    datamean = mean(data, 1); 
    datasd = std(data, [], 1); 
    zscored_data = (data-repmat(datamean, size(data,1), 1))./repmat(datasd, size(data,1), 1); 
    
    % Correlation matrix
    figure 
    imagesc(corr(zscored_data, 'Type', 'Pearson'))
    axis square
    xticks(1:size(data,2))
    xticklabels(var_names)
    xtickangle(60)
    yticks(1:size(data,2))
    yticklabels(var_names)
    title('Correlation matrix - z-scored data')
    c = colorbar; c.Label.String = 'Correlation coefficient';
    
    % Run PCA
    [coeff,score,latent,tsquared,explained,mu] = pca(zscored_data,'rows','pairwise');
    
    % Scree plot 
    figure 
    subplot(2,2,[1 2]), hold on
    plot(1:(length(latent)), latent, 'ko-')
    plot(get(gca, 'xlim'), [1 1]*mean(latent), 'k--')
    xlabel('Principal components')
    ylabel('Eigen value')
    title('Scree plot'), hold off 
    
    subplot(2, 2, [3, 4])
    plot(1:(length(latent)), cumsum(explained), 'ko-')
    xlabel('Principal components')
    ylabel('Cumulative explained variance (%)')
    
    % Quality of representation  (cos2)
    loadings = coeff .* repmat(sqrt(latent)', size(coeff,1), 1); 
    cos2 = loadings.^2; 
    
    figure
    subplot(1,2,1)
    imagesc(cos2)
    axis square, colormap jet
    title('Quality of representation (cos^2)')
    xticks(1:size(coeff, 2))
    xlabel('Principal components')
    yticks(1:size(coeff, 2))
    yticklabels(var_names)
    caxis([0 1])
    c = colorbar; c.Label.String = 'cos^2';
    
    % Loadings 
    subplot(1,2,2) 
    imagesc(loadings)
    axis square, colormap jet
    title('Loadings')
    xticks(1:size(coeff, 2))
    xlabel('Principal components')
    yticks(1:size(coeff, 2))
    yticklabels(var_names)
    caxis([-1 1])
    c = colorbar; c.Label.String = 'Corr. original observations - PC scores';

    % Correlation circle 
    th=0:pi/100:2*pi; xcircle= cos(th); ycircle = sin(th); 
    for pci = 1:n2plot
        for pcj = pci+1:n2plot
            figure
            plot(xcircle,ycircle,'k--'); hold on 
            plot([0 0], [-1 1], 'k-'); plot([-1 1], [0 0], 'k-') 
            scatter(loadings(:, pci), loadings(:,pcj), 'bo')
            plot([zeros(size(loadings, 1), 1), loadings(:, pci)]', ...
                [zeros(size(loadings, 1), 1), loadings(:,pcj)]', 'b:')
            text(loadings(:, pci)+.05, loadings(:,pcj), var_names, 'FontSize',9);
            title(['PC' num2str(pci) ' vs PC' num2str(pcj)])
            xlabel(['PC' num2str(pci)]); ylabel(['PC' num2str(pcj)])
            set(gca, 'xlim', [-1 1]); set(gca, 'ylim', [-1 1])
            xticks(-1:.2:1); yticks(-1:.2:1); 
            grid on, axis square, hold off 
        end
    end

    % Contributions
    contrib_pct = cos2./repmat(sum(cos2,1), size(cos2, 1), 1)*100; 
    figure
    plot_count = 1; 
    for pci = 1:n2plot
        subplot(ceil(sqrt(n2plot)), ceil(sqrt(n2plot)), plot_count)
        plot_count = plot_count+1; 
        pc_contrib = [contrib_pct(:, pci), (1:size(contrib_pct, 1))']; 
        pc_contrib = sortrows(pc_contrib, 'descend');
        bar(pc_contrib(:,1), 'EdgeColor', 'none')
        xticks(1:size(contrib_pct, 1))
        xticklabels(var_names(pc_contrib(:,2)))
        xtickangle(60)
        ylabel('Contrib. to PC (%)')
        title(['Contributions to PC ' num2str(pci)])
    end
    
end 