function TDepPlot(struct1, struct2, metricType, titleStr, binaural)
    % Plot the time-dependent overall sound quality comparison
    
    switch metricType
        case 'tonality'
            cmap = load('cmap_plasma.txt');
            unit = "Tonality, tu_{HMS}";
            metric1 = struct1.TonalTDep;
            metric2 = struct2.tonalityTDep;
        case 'loudness'
            cmap = load('cmap_viridis.txt');
            unit = "Loudness, sone_{HMS}";
            if binaural
                metric1 = struct1.LoudTDepBin;
                metric2 = struct2.loudnessTDepBin;
            else
                metric1 = struct1.LoudTDep;
                metric2 = struct2.loudnessTDep;
            end
        case 'roughness'
            cmap = load('cmap_inferno.txt');
            unit = "Roughness, asper_{HMS}";
            if binaural
                metric1 = struct1.RoughTDepBin;
                metric2 = struct2.roughnessTDepBin;
            else
                metric1 = struct1.RoughTDep;
                metric2 = struct2.roughnessTDep;
            end
    end

    if size(metric2, 2) == 2
        fig = figure('Position', [200, 200, 900, 300]);
        tiledlayout(fig, 1, 2);
        titleStr = [strcat(titleStr, " left"), strcat(titleStr, " right")];
    elseif size(metric2, 2) == 1
        fig = figure('Position', [200, 200, 450, 300]);
        tiledlayout(fig, 1, 1);
        if binaural && ~(strcmp(metricType, 'tonality'))
            titleStr = strcat(titleStr, " binaural");
        end
    else
        error("Check your inputs!")
    end

    movegui(fig, 'center');

    cmap1 = 166;
    cmap2 = 34;

    for ii = 1:(size(metric2, 2))
    
        ax = nexttile(ii);
        plot(metric1(:, 1), metric1(:, ii + 1),...
             'Color',  cmap(cmap1, :),...
             'LineWidth', 1, 'LineStyle', '-',...
             'DisplayName', "ArtemiS");
        hold on
        plot(struct2.timeOut, metric2(:, ii),...
             'Color',  cmap(cmap2, :),...
             'LineWidth', 1.5, 'LineStyle', ':',...
             'DisplayName', "refmap");
        hold off
        ax.XLim = [struct2.timeOut(1), struct2.timeOut(end)...
                   + (struct2.timeOut(2) - struct2.timeOut(1))];
        ax.XLabel.String = "Time, s";
        ax.YLim = [0, 1.1*ceil(max(metric2(:, ii))*10)/10];
        ax.YLabel.String = unit;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridAlpha = 0.15;
        ax.FontName = 'Arial';
        ax.FontSize = 10;
        ax.Title.String = titleStr(ii);
        ax.Title.FontWeight = 'normal';
        ax.Title.FontSize = 11;
        lgd = legend('Location', 'best', 'FontSize', 8);
        lgd.FontSize = 10;
    end
end