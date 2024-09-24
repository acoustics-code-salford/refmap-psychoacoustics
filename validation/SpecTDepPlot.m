function SpecTDepPlot(struct1, struct2, metricType, titleStr, binaural)
    % Plot the time-dependent specific sound quality comparison

    switch metricType
        case 'tonality'
            cmap = load('cmap_plasma.txt');
            unit = {"Specific tonality,"; "tu_{HMS}/Bark_{HMS}"};
            metric1 = struct1.TonalSpecTDep;
            metric2 = struct2.specTonality;
        case 'loudness'
            cmap = load('cmap_viridis.txt');
            unit = {"Specific loudness,"; "sone_{HMS}/Bark_{HMS}"};
            if binaural
                metric1 = struct1.LoudSpecTDepBin;
                metric2 = struct2.specLoudnessBin;
            else
                metric1 = struct1.LoudSpecTDep;
                metric2 = struct2.specLoudness;
            end
        case 'roughness'
            cmap = load('cmap_inferno.txt');
            unit = {"Specific roughness,"; "asper_{HMS}/Bark_{HMS}"};
            if binaural
                metric1 = struct1.RoughSpecTDepBin;
                metric2 = struct2.specRoughnessBin;
            else
                metric1 = struct1.RoughSpecTDep;
                metric2 = struct2.specRoughness;
            end
    end

    if size(metric2, 3) == 2
        fig = figure('Position', [200, 200, 1000, 550]);
        tl = tiledlayout(fig, 2, 2);
        axTitleStr = ["ArtemiS left", "refmap left", "ArtemiS right", "refmap right"];
        metric1 = reshape([metric1(1:size(metric1, 1)/2, :), metric1(size(metric1, 1)/2 + 1:end, :)], [size(metric1, 1)/2, size(metric1, 2), 2]);
        axs = [1, 3, 2, 4];
    elseif size(metric2, 3) == 1
        fig = figure('Position', [200, 200, 1000, 300]);
        tl = tiledlayout(fig, 1, 2);
        axTitleStr = ["ArtemiS", "refmap"];
        if binaural && ~(strcmp(metricType, 'tonality'))
            titleStr = strcat(titleStr, " binaural");
        end
        axs = [1, 2];
    else
        error("Check your inputs!")
    end

    movegui(fig, 'center');

    for ii = 1:length(axs)
        ax = nexttile(axs(ii));
        
        if ii == 1 || ii == 3
            if ii == 1
                surf(ax, metric1(2:end, 1, 1), struct2.bandCentreFreqs,...
                     permute(metric1(2:end, 2:end, 1), [2, 1, 3]),...
                     'EdgeColor', 'none', 'FaceColor', 'interp',...
                     'DisplayName', "ArtemiS");
            else
                surf(ax, metric1(2:end, 1, 2), struct2.bandCentreFreqs,...
                     permute(metric1(2:end, 2:end, 2), [2, 1, 3]),...
                     'EdgeColor', 'none', 'FaceColor', 'interp',...
                     'DisplayName', "ArtemiS");
            end

        else
            if ii == 2
                surf(ax, struct2.timeOut, struct2.bandCentreFreqs,...
                     permute(metric2(:, :, 1), [2, 1, 3]),...
                     'EdgeColor', 'none', 'FaceColor', 'interp',...
                     'DisplayName', "refmap");
            else
                surf(ax, struct2.timeOut, struct2.bandCentreFreqs,...
                     permute(metric2(:, :, 2), [2, 1, 3]),...
                     'EdgeColor', 'none', 'FaceColor', 'interp',...
                     'DisplayName', "refmap");
            end
        end
        
        view(2);

        ax.XLim = [struct2.timeOut(1), struct2.timeOut(end) + (struct2.timeOut(2) - struct2.timeOut(1))];
        ax.YLim = [struct2.bandCentreFreqs(1), struct2.bandCentreFreqs(end)];
        ax.CLim = [0, ceil(max(metric2, [], 'all')*10)/10];
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
        colormap(cmap);
        h = colorbar;
        set(get(h,'label'),'string', unit);
    end

    tl.Title.String = titleStr;
    tl.Title.FontSize = 11;
end