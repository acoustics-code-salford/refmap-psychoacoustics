function PyDiffSpecTDepPlot(pyMetric, mlabMetric, metricType, titleStr,...
                            figPath, figName)
    % Plot the time-dependent specific sound quality differences between
    % Python and MATLAB implementations

    switch metricType
        case 'tonality'
            cmap = load('cmap_plasma.txt');
            unit = {"Specific tonality (difference),"; "\Delta tu_{SHM}/Bark_{SHM}"};
            metricPy = pyMetric.specTonality;
            metricMlab = mlabMetric.specTonality;
            
        case 'loudness'
            cmap = load('cmap_viridis.txt');
            unit = {"Specific loudness (difference),"; "\Delta sone_{SHM}/Bark_{SHM}"};
            metricPy = pyMetric.specLoudness;
            metricMlab = mlabMetric.specLoudness;
        case 'roughness'
            cmap = load('cmap_inferno.txt');
            unit = {"Specific roughness (difference),"; "\Delta asper_{SHM}/Bark_{SHM}"};
            metricPy = pyMetric.specRoughness;
            metricMlab = mlabMetric.specRoughness;
    end

    metricDiff = metricPy - metricMlab;

    fig = figure('Position', [200, 200, 500, 700]);
    tl = tiledlayout(fig, 2, 1);
    axTitleStr = ["Left channel", "Right channel"];

    movegui(fig, 'center');

    for ii = 1:2
        ax = nexttile(ii);
        surf(ax, pyMetric.timeOut, pyMetric.bandCentreFreqs,...
             permute(metricDiff(:, :, ii), [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');        
        view(2);

        ax.XLim = [pyMetric.timeOut(1), pyMetric.timeOut(end)...
                   + (pyMetric.timeOut(2) - pyMetric.timeOut(1))];
        ax.YLim = [pyMetric.bandCentreFreqs(1), pyMetric.bandCentreFreqs(end)];
        % ax.CLim = [-0.005, 0.005];
        ax.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
        ax.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                         "8k", "16k"]; 
        ax.YScale = 'log';
        ax.YLabel.String = "Frequency, Hz";
        ax.XLabel.String = "Time, s";
        ax.FontName = 'Arial';
        ax.FontSize = 10;
        ax.Title.String = axTitleStr(ii);
        ax.Title.FontWeight = 'normal';
        ax.Title.FontSize = 10;
        ax.Toolbar.Visible = 'off';
        colormap(cmap);
        h = colorbar;
        set(get(h,'label'),'string', unit);
    end

    tl.Title.String = titleStr;
    tl.Title.FontSize = 11;

    if isstring(figPath)
        exportgraphics(fig, fullfile(figPath, figName + ".png"), 'Resolution', 300)
    end
end