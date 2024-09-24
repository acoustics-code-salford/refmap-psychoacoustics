function SpecTAggPlot(struct1, struct2, metricType, titleStr, binaural)
    % Plot the time-aggregated specific sound quality comparison
    
    switch metricType
        case 'tonality'
            cmap = load('cmap_plasma.txt');
            unit = "Specific tonality, tu_{HMS}/Bark_{HMS}";
            metric1 = struct1.TonalSpec;
            metric2 = struct2.specTonalityAvg;
        case 'loudness'
            cmap = load('cmap_viridis.txt');
            unit = "Specific loudness, sone_{HMS}/Bark_{HMS}";
            if binaural
                metric1 = struct1.LoudSpecBin;
                metric2 = struct2.specLoudnessPowAvgBin;
            else
                metric1 = struct1.LoudSpec;
                metric2 = struct2.specLoudnessPowAvg;
            end
        case 'roughness'
            cmap = load('cmap_inferno.txt');
            unit = "Specific roughness, asper_{HMS}/Bark_{HMS}";
            if binaural
                metric1 = struct1.RoughSpecBin;
                metric2 = struct2.specRoughnessAvgBin;
            else
                metric1 = struct1.RoughSpec;
                metric2 = struct2.specRoughnessAvg;
            end
    end
    
    if size(metric2, 2) == 2
        fig = figure('Position', [200, 200, 1000, 325]);
        tiledlayout(fig, 1, 2);
        titleStr = [strcat(titleStr, " left"), strcat(titleStr, " right")];
    elseif size(metric2, 2) == 1
        fig = figure('Position', [200, 200, 500, 325]);
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
        bar(linspace(0.5, 26.5, 53),...
            metric1(:, ii + 1),...
            'EdgeColor',  cmap(cmap1, :),...
            'FaceColor', 'none',...
            'LineWidth', 1, 'LineStyle', '-',...
            'DisplayName', "ArtemiS");
        hold on
        bar(linspace(0.5, 26.5, 53),...
            metric2(:, ii),...
            'EdgeColor',  cmap(cmap2, :),...
            'FaceColor',  'none',...
            'LineWidth', 0.5, 'LineStyle', '--',...
            'DisplayName', "refmap");
        hold off
        ax.XTick = linspace(0.5, 26.5, 27);
        ax.XLabel.String = "Critical band rate, Bark_{HMS}";
        ax.YLim = [0, 1.1*ceil(max(metric2(:, ii))*10)/10];
        ax.YLabel.String = unit;
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
