function acousticHMSValidation
% acousticHMSValidation
%
% Compares ECMA-418-2:2022 implementation with values obtained using
% commercially-available software with reference signals.
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 19/08/2024
% Date last modified: 24/09/2024
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
% Checked by:
% Date last checked:
%
% input parameters
%
%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))

%% Import reference data
refpath = fullfile("validation", "ECMA-418-2_2022", "reference");

% Loudness

sine_1kHz_40dB.LoudTDep = readmatrix(fullfile(refpath, "sine_1kHz_40dB.Loud_HMS_TDep_LR.asc"), 'FileType', 'text');
sine_1kHz_40dB.LoudSpec = readmatrix(fullfile(refpath, "sine_1kHz_40dB.LoudSpec_HMS_LR.asc"), 'FileType', 'text');
sine_1kHz_40dB.LoudSpecTDep = readmatrix(fullfile(refpath, "sine_1kHz_40dB.LoudSpec_HMS_TDep_LR.asc"), 'FileType', 'text');

sine_1kHz_70Hz_60dB.LoudTDep = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.Loud_HMS_TDep_LR.asc"), 'FileType', 'text');
sine_1kHz_70Hz_60dB.LoudSpec = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.LoudSpec_HMS_LR.asc"), 'FileType', 'text');
sine_1kHz_70Hz_60dB.LoudSpecTDep = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.LoudSpec_HMS_TDep_LR.asc"), 'FileType', 'text');

BusyStreet1_0530_0600.LoudTDep = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.Loud_HMS_TDep_LR.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.LoudSpec = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.LoudSpec_HMS_LR.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.LoudSpecTDep = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.LoudSpec_HMS_TDep_LR.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.LoudTDepBin = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.Loud_HMS_TDep_Bin.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.LoudSpecBin = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.LoudSpec_HMS_Bin.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.LoudSpecTDepBin = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.LoudSpec_HMS_TDep_Bin.asc"), 'FileType', 'text');

% Tonality

sine_1kHz_40dB.TonalTDep = readmatrix(fullfile(refpath, "sine_1kHz_40dB.Tonal_HMS_TDep_LR.asc"), 'FileType', 'text');
sine_1kHz_40dB.TonalSpec = readmatrix(fullfile(refpath, "sine_1kHz_40dB.TonalSpec_HMS_LR.asc"), 'FileType', 'text');
sine_1kHz_40dB.TonalSpecTDep = readmatrix(fullfile(refpath, "sine_1kHz_40dB.TonalSpec_HMS_TDep_LR.asc"), 'FileType', 'text');

sine_1kHz_70Hz_60dB.TonalTDep = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.Tonal_HMS_TDep_LR.asc"), 'FileType', 'text');
sine_1kHz_70Hz_60dB.TonalSpec = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.TonalSpec_HMS_LR.asc"), 'FileType', 'text');
sine_1kHz_70Hz_60dB.TonalSpecTDep = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.TonalSpec_HMS_TDep_LR.asc"), 'FileType', 'text');

BusyStreet1_0530_0600.TonalTDep = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.Tonal_HMS_TDep_LR.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.TonalSpec = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.TonalSpec_HMS_LR.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.TonalSpecTDep = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.TonalSpec_HMS_TDep_LR.asc"), 'FileType', 'text');

% Roughness

sine_1kHz_70Hz_60dB.RoughTDep = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.Rough_HMS_TDep_LR.asc"), 'FileType', 'text');
sine_1kHz_70Hz_60dB.RoughSpec = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.RoughSpec_HMS_LR.asc"), 'FileType', 'text');
sine_1kHz_70Hz_60dB.RoughSpecTDep = readmatrix(fullfile(refpath, "sine_1kHz_70Hz_60dB.RoughSpec_HMS_TDep_LR.asc"), 'FileType', 'text');

BusyStreet1_0530_0600.RoughTDep = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.Rough_HMS_TDep_LR.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.RoughSpec = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.RoughSpec_HMS_LR.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.RoughSpecTDep = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.RoughSpec_HMS_TDep_LR.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.RoughTDepBin = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.Rough_HMS_TDep_Bin.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.RoughSpecBin = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.RoughSpec_HMS_Bin.asc"), 'FileType', 'text');
BusyStreet1_0530_0600.RoughSpecTDepBin = readmatrix(fullfile(refpath, "BusyStreet1_0530-0600.RoughSpec_HMS_TDep_Bin.asc"), 'FileType', 'text');

% single values
validation_ECMA_418_2_2_2022_LR = readcell(fullfile(refpath, "validation_ECMA-418-2_2022_LR.xlsx"),...
                                            'NumHeaderLines', 0);

sine_1kHz_40dB.Loudness = validation_ECMA_418_2_2_2022_LR{2, 3};
sine_1kHz_70Hz_60dB.Loudness = validation_ECMA_418_2_2_2022_LR{3, 3};
BusyStreet1_0530_0600.Loudness = [validation_ECMA_418_2_2_2022_LR{4, 3},...
                                  validation_ECMA_418_2_2_2022_LR{5, 3}];

sine_1kHz_40dB.Tonality = validation_ECMA_418_2_2_2022_LR{2, 5};
sine_1kHz_70Hz_60dB.Tonality = validation_ECMA_418_2_2_2022_LR{3, 5};
BusyStreet1_0530_0600.Tonality = [validation_ECMA_418_2_2_2022_LR{4, 5},...
                                  validation_ECMA_418_2_2_2022_LR{5, 5}];

sine_1kHz_40dB.Roughness = validation_ECMA_418_2_2_2022_LR{2, 4};
sine_1kHz_70Hz_60dB.Roughness = validation_ECMA_418_2_2_2022_LR{3, 4};
BusyStreet1_0530_0600.Roughness = [validation_ECMA_418_2_2_2022_LR{4, 4},...
                                   validation_ECMA_418_2_2_2022_LR{5, 4}];

validation_ECMA_418_2_2_2022_Bin = readcell(fullfile(refpath, "validation_ECMA-418-2_2022_Bin.xlsx"),...
                                            'NumHeaderLines', 0);

BusyStreet1_0530_0600.LoudnessBin = validation_ECMA_418_2_2_2022_Bin{2, 2};
BusyStreet1_0530_0600.RoughnessBin = validation_ECMA_418_2_2_2022_Bin{2, 3};

% concatenate results

signalLabs = string({validation_ECMA_418_2_2_2022_LR{:, 1}}) + " " + string({validation_ECMA_418_2_2_2022_LR{:, 2}});
signalLabs = signalLabs(2:end);
signalLabs = eraseBetween(signalLabs, " ", "Ch");
signalLabs = [signalLabs, extractBefore(signalLabs(end), " ") + " Bin"];

loudSingles = horzcat([sine_1kHz_40dB.Loudness, sine_1kHz_70Hz_60dB.Loudness],...
                      BusyStreet1_0530_0600.Loudness, BusyStreet1_0530_0600.LoudnessBin);

tonalSingles = horzcat([sine_1kHz_40dB.Tonality, sine_1kHz_70Hz_60dB.Tonality],...
                       BusyStreet1_0530_0600.Tonality);

roughSingles = horzcat([sine_1kHz_40dB.Roughness, sine_1kHz_70Hz_60dB.Roughness],...
               BusyStreet1_0530_0600.Roughness, BusyStreet1_0530_0600.RoughnessBin);

%% Import audio data

[signal1, fs1] = audioread("sine_1kHz_40dB.wav");
[signal2, fs2] = audioread("sine_1kHz_70Hz_60dB.wav");
[signal3, fs3] = audioread("BusyStreet1_0530-0600.wav");

%% Calculate sound qualities

tonalityHMS1 = acousticHMSTonality(signal1, fs1, 1, false);
tonalityHMS2 = acousticHMSTonality(signal2, fs2, 1, false);
tonalityHMS3 = acousticHMSTonality(signal3, fs3, 1, false);

loudnessHMS1 = acousticHMSLoudnessFromComponent(tonalityHMS1.specTonalLoudness,...
                                                tonalityHMS1.specNoiseLoudness,...
                                                false, false);
loudnessHMS1full = acousticHMSLoudness(signal1, fs1, 1, false, false);
loudnessHMS2 = acousticHMSLoudnessFromComponent(tonalityHMS2.specTonalLoudness,...
                                                tonalityHMS2.specNoiseLoudness,...
                                                false, false);
loudnessHMS3 = acousticHMSLoudnessFromComponent(tonalityHMS3.specTonalLoudness,...
                                                tonalityHMS3.specNoiseLoudness,...
                                                false, true);

roughnessHMS1 = acousticHMSRoughness(signal1, fs1, 1, false, false);
roughnessHMS2 = acousticHMSRoughness(signal2, fs2, 1, false, false);
roughnessHMS3 = acousticHMSRoughness(signal3, fs3, 1, false, true);

loudSinglesAll = vertcat(loudSingles, horzcat([loudnessHMS1.loudnessPowAvg,...
                                               loudnessHMS2.loudnessPowAvg],...
                                              loudnessHMS3.loudnessPowAvg,...
                                              loudnessHMS3.loudnessPowAvgBin));

tonalSinglesAll = vertcat(tonalSingles, horzcat([tonalityHMS1.tonalityAvg,...
                                                 tonalityHMS2.tonalityAvg],...
                                                tonalityHMS3.tonalityAvg));

roughSinglesAll = vertcat(roughSingles, horzcat([roughnessHMS1.roughness90Pc,...
                                                 roughnessHMS2.roughness90Pc],...
                                                roughnessHMS3.roughness90Pc,...
                                                roughnessHMS3.roughness90PcBin));


%% Results comparison plots

% path for saving figures
figpath = fullfile("validation", "ECMA-418-2_2022", "results");

% figure format
figformat = 'png';

% Tonality
% --------

% sine_1kHz_40dB.wav
% Time-dependent
TDepPlot(sine_1kHz_40dB, tonalityHMS1, 'tonality', "sine\_1kHz\_40dB.wav",...
         false, figpath, "tonalHMSTDepSine1kHz40dB", figformat);
SpecTDepPlot(sine_1kHz_40dB, tonalityHMS1, 'tonality', "sine\_1kHz\_40dB.wav",...
         false, figpath, "tonalHMSSpecTDepSine1kHz40dB");

% Time-aggregated
SpecTAggPlot(sine_1kHz_40dB, tonalityHMS1, 'tonality', "sine\_1kHz\_40dB.wav", false, NaN, NaN);

% sine_1kHz_70Hz_60dB.wav
% Time-dependent
TDepPlot(sine_1kHz_70Hz_60dB, tonalityHMS2, 'tonality', "sine\_1kHz\_70Hz\_60dB.wav", false, NaN, NaN);
SpecTDepPlot(sine_1kHz_70Hz_60dB, tonalityHMS2, 'tonality', "sine\_1kHz\_70Hz\_60dB.wav", false, NaN, NaN);

% Time-aggregated
SpecTAggPlot(sine_1kHz_70Hz_60dB, tonalityHMS2, 'tonality', "sine\_1kHz\_70Hz\_60dB.wav", false, NaN, NaN);

% BusyStreet1_0530-0600.wav
% Time-dependent
TDepPlot(BusyStreet1_0530_0600, tonalityHMS3, 'tonality', "BusyStreet1\_0530-0600.wav",...
         false, figpath, "tonalHMSTDepBusySt", figformat);
SpecTDepPlot(BusyStreet1_0530_0600, tonalityHMS3, 'tonality', "BusyStreet1\_0530-0600.wav",...
             false, figpath, "tonalHMSSpecTDepBusySt");

% Time-aggregated
SpecTAggPlot(BusyStreet1_0530_0600, tonalityHMS3, 'tonality', "BusyStreet1\_0530-0600.wav",...
             false, figpath, "tonalHMSSpecTAggBusySt", figformat);

% Single values
fg = figure('Position', [200, 200, 450, 300]);
movegui(fg, 'center');
br = bar(strrep(signalLabs(1:end - 1), "_", " "), tonalSinglesAll);
colororder('dye')
ylabel("Tonality, tu_{HMS}")
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'GridAlpha', 0.15);
legend(["ArtemiS", "refmap"], 'Location', 'northeast')
for bb = 1:length(br)
    text(br(bb).XEndPoints, br(bb).YEndPoints, string(round(br(bb).YData, 2)),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', 8)
end
exportgraphics(fg, fullfile(figpath, "tonalHMSsingles.pdf"), 'ContentType', 'vector')
exportgraphics(fg, fullfile(figpath, "tonalHMSsingles.png"), 'Resolution', 300)

% Loudness
% --------

% sine_1kHz_40dB.wav
% Time-dependent (from component)
TDepPlot(sine_1kHz_40dB, loudnessHMS1, 'loudness', "sine\_1kHz\_40dB.wav (from component)",...
         false, figpath, "loudHMSTDepSine1kHz40dB", figformat);
SpecTDepPlot(sine_1kHz_40dB, loudnessHMS1, 'loudness', "sine\_1kHz\_40dB.wav (from component)",...
             false, figpath, "loudHMSSpecTDepSine1kHz40dB");

% Time-aggregated (from component)
SpecTAggPlot(sine_1kHz_40dB, loudnessHMS1, 'loudness',...
             "sine\_1kHz\_40dB.wav (from component)", false, NaN, NaN);

% Time-dependent (full function)
TDepPlot(sine_1kHz_40dB, loudnessHMS1full, 'loudness',...
         "sine\_1kHz\_40dB.wav (full function)", false, NaN, NaN);
SpecTDepPlot(sine_1kHz_40dB, loudnessHMS1full, 'loudness',...
             "sine\_1kHz\_40dB.wav (full function)", false, NaN, NaN);

% Time-aggregated (full function)
SpecTAggPlot(sine_1kHz_40dB, loudnessHMS1full, 'loudness',...
             "sine\_1kHz\_40dB.wav (full function)", false, NaN, NaN);

% sine_1kHz_70Hz_60dB.wav
% Time-dependent
TDepPlot(sine_1kHz_70Hz_60dB, loudnessHMS2, 'loudness',...
         "sine\_1kHz\_70Hz\_60dB.wav", false, NaN, NaN);
SpecTDepPlot(sine_1kHz_70Hz_60dB, loudnessHMS2, 'loudness',...
             "sine\_1kHz\_70Hz\_60dB.wav", false, NaN, NaN);

% Time-aggregated
SpecTAggPlot(sine_1kHz_70Hz_60dB, loudnessHMS2, 'loudness',...
             "sine\_1kHz\_70Hz\_60dB.wav", false, NaN, NaN);

% BusyStreet1_0530-0600.wav
% Time-dependent (separate)
TDepPlot(BusyStreet1_0530_0600, loudnessHMS3, 'loudness', "BusyStreet1\_0530-0600.wav",...
         false, figpath, "loudHMSTDepBusySt", figformat);
SpecTDepPlot(BusyStreet1_0530_0600, loudnessHMS3, 'loudness', "BusyStreet1\_0530-0600.wav",...
             false, figpath, "loudHMSSpecTDepBusySt");

% Time-aggregated (separate)
SpecTAggPlot(BusyStreet1_0530_0600, loudnessHMS3, 'loudness', "BusyStreet1\_0530-0600.wav",...
             false, figpath, "loudHMSSpecTAggBusySt", figformat);

% Time-dependent (binaural)
TDepPlot(BusyStreet1_0530_0600, loudnessHMS3, 'loudness',...
         "BusyStreet1\_0530-0600.wav", true, NaN, NaN);
SpecTDepPlot(BusyStreet1_0530_0600, loudnessHMS3, 'loudness',...
             "BusyStreet1\_0530-0600.wav", true, figpath, "loudHMSSpecTDepBinBusySt");

% Time-aggregated (binaural)
SpecTAggPlot(BusyStreet1_0530_0600, loudnessHMS3, 'loudness',...
             "BusyStreet1\_0530-0600.wav", true, NaN, NaN);

% Single values
fg = figure('Position', [200, 200, 500, 300]);
movegui(fg, 'center');
br = bar(strrep(signalLabs, "_", " "), loudSinglesAll);
colororder('earth')
ylabel("Loudness, sone_{HMS}")
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'GridAlpha', 0.15);
legend(["ArtemiS", "refmap"], 'Location', 'northwest')
for bb = 1:length(br)
    text(br(bb).XEndPoints, br(bb).YEndPoints, string(round(br(bb).YData, 2)),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', 7)
end
exportgraphics(fg, fullfile(figpath, "loudHMSsingles.pdf"), 'ContentType', 'vector')
exportgraphics(fg, fullfile(figpath, "loudHMSsingles.png"), 'Resolution', 300)

% Roughness
% ---------

% sine_1kHz_70Hz_60dB.wav
% Time-dependent
TDepPlot(sine_1kHz_70Hz_60dB, roughnessHMS2, 'roughness', "sine\_1kHz\_70Hz\_60dB.wav",...
         false, figpath, "roughHMSTDepSine1kHz70Hz60dB", figformat);
SpecTDepPlot(sine_1kHz_70Hz_60dB, roughnessHMS2, 'roughness', "sine\_1kHz\_70Hz\_60dB.wav",...
             false, figpath, "roughHMSSpecTDepSine1kHz70Hz60dB");

% Time-aggregated
SpecTAggPlot(sine_1kHz_70Hz_60dB, roughnessHMS2, 'roughness',...
             "sine\_1kHz\_70Hz\_60dB.wav", false, NaN, NaN);

% BusyStreet1_0530-0600.wav
% Time-dependent (separate)
TDepPlot(BusyStreet1_0530_0600, roughnessHMS3, 'roughness', "BusyStreet1\_0530-0600.wav",...
         false, figpath, "roughHMSTDepBusySt", figformat);
SpecTDepPlot(BusyStreet1_0530_0600, roughnessHMS3, 'roughness', "BusyStreet1\_0530-0600.wav",...
             false, figpath, "roughHMSSpecTDepBusySt");

% Time-aggregated (separate)
SpecTAggPlot(BusyStreet1_0530_0600, roughnessHMS3, 'roughness', "BusyStreet1\_0530-0600.wav",...
             false, figpath, "roughHMSSpecTAggBusySt", figformat);

% Time-dependent (binaural)
TDepPlot(BusyStreet1_0530_0600, roughnessHMS3, 'roughness',...
         "BusyStreet1\_0530-0600.wav", true, NaN, NaN);
SpecTDepPlot(BusyStreet1_0530_0600, roughnessHMS3, 'roughness', "BusyStreet1\_0530-0600.wav",...
             true, figpath, "roughHMSSpecTDepBinBusySt");

% Time-aggregated (binaural)
SpecTAggPlot(BusyStreet1_0530_0600, roughnessHMS3, 'roughness',...
             "BusyStreet1\_0530-0600.wav", true, NaN, NaN);

% Single values
fg = figure('Position', [200, 200, 500, 300]);
movegui(fg, 'center');
br = bar(strrep(signalLabs, "_", " "), roughSinglesAll);
colororder('gem')
ylabel("Roughness, asper_{HMS}")
set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'GridAlpha', 0.15);
legend(["ArtemiS", "refmap"], 'Location', 'northeast')
for bb = 1:length(br)
    text(br(bb).XEndPoints, br(bb).YEndPoints, string(round(br(bb).YData, 2)),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', 7)
end
exportgraphics(fg, fullfile(figpath, "roughHMSsingles.pdf"), 'ContentType', 'vector')
exportgraphics(fg, fullfile(figpath, "roughHMSsingles.png"), 'Resolution', '300')

%% Plotting functions

function TDepPlot(struct1, struct2, metricType, titleStr, binaural,...
                  figPath, figName, figFormat)
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
        ax.GridLineStyle = '--';
        ax.GridAlpha = 0.15;
        ax.FontName = 'Arial';
        ax.FontSize = 10;
        ax.Title.String = titleStr(ii);
        ax.Title.FontWeight = 'normal';
        ax.Title.FontSize = 11;
        lgd = legend('Location', 'best', 'FontSize', 8);
        lgd.FontSize = 10;
    end

    if isstring(figPath)
        if strcmpi(figFormat, "pdf")
            exportgraphics(fig, fullfile(figPath, figName + ".pdf"), 'ContentType', 'vector')
        elseif strcmpi(figFormat, "png")
            exportgraphics(fig, fullfile(figPath, figName + ".png"), 'Resolution', 300)
        end
    end
end


function SpecTAggPlot(struct1, struct2, metricType, titleStr, binaural,...
                      figPath, figName, figFormat)
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
        ax.GridLineStyle = '--';
        ax.GridAlpha = 0.15;
        ax.FontName = 'Arial';
        ax.FontSize = 10;
        ax.Title.String = titleStr(ii);
        ax.Title.FontWeight = 'normal';
        ax.Title.FontSize = 11;
        lgd = legend('Location', 'best', 'FontSize', 8);
        lgd.FontSize = 10;
    end

    if isstring(figPath)
        if strcmpi(figFormat, "pdf")
            exportgraphics(fig, fullfile(figPath, figName + ".pdf"), 'ContentType', 'vector')
        elseif strcmpi(figFormat, "png")
            exportgraphics(fig, fullfile(figPath, figName + ".png"), 'Resolution', 300)
        end
    end
end

function SpecTDepPlot(struct1, struct2, metricType, titleStr, binaural,...
                      figPath, figName)
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

    if isstring(figPath)
        exportgraphics(fig, fullfile(figPath, figName + ".png"), 'Resolution', 300)
    end
end

end