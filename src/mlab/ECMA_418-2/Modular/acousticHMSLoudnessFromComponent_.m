function loudnessHMS = acousticHMSLoudnessFromComponent_(specTonalLoudness,...
                                                         specNoiseLoudness,...
                                                         outplot, ecma, binaural)
% [loudnessPowAvg, loudnessTimeVar, specificLoudness, 
%  specificLoudnessPowAvg, bandCentreFreqs]
%  = acousticHMSLoudnessFromComponent_(specificTonalLoudness,
%                                     specificNoiseLoudness,
%                                     outplot, ecma, binaural)
%
% Returns loudness values according to ECMA-418-2:2022 (using the Hearing
% Model of Sottek) for input component specific tonal loudness and specific
% noise loudness, obtained using acousticHMSTonality_().
%
% Inputs
% ------
% specTonalLoudness : 2D or 3D matrix
%                     the specific tonal loudness values calculated for
%                     a sound pressure signal (single mono or single
%                     stereo audio)
%
% specNoiseLoudness : 2D or 3D matrix
%                     the specific noise loudness values calculated for
%                     a sound pressure signal (single mono or single
%                     stereo audio)
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% ecma : Boolean true/false (default: true)
%        flag indicating whether to maintain strict standard adherence to
%        ECMA-418-2:2022 Equation 40, or otherwise to use an alternative
%        that provides closer time-alignment of the time-dependent loudness
%        with the original signal
%
% binaural : Boolean true/false (default: true)
%            flag indicating whether to output binaural loudness for stereo
%            input signal.
% 
% Returns
% -------
% For each 'channel' in the input matrices:
%
% loudnessPowAvg : number or vector
%                  (power)-average (overall) loudness value
% 
% loudnessTimeVar : vector or 2D matrix
%                   time-dependent overall loudness values
%
% specificLoudness : 2D or 3D matrix
%                    time-dependent specific loudness values in each
%                    half-critical band rate scale width
%
% specificLoudnessPowAvg : vector or 2D matrix
%                          time-(power)-averaged specific loudness values
%                          in each half-critical band rate scale width
%
% bandCentreFreqs : vector
%                   centre frequencies corresponding with each half-Bark
%                   critical band rate scale width
%
% Assumptions
% -----------
% The input matrices are ECMA-418-2:2022 specific tonal and specific noise
% loudness, with dimensions orientated as [half-Bark bands, time blocks,
% signal channels]
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 22/08/2023
% Date last modified: 10/07/2024
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        specTonalLoudness (:, :, :) double {mustBeReal}
        specNoiseLoudness (:, :, :) double {mustBeReal}
        outplot {mustBeNumericOrLogical} = false
        ecma {mustBeNumericOrLogical} = true
        binaural {mustBeNumericOrLogical} = true
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))

%% Input checks
% Check the size of the input matrices (must match)
if ~isequal(size(specTonalLoudness), size(specNoiseLoudness))
    error('Error: Input loudness matrix sizes must match')
end

% Check the channel number of the input data
if size(specTonalLoudness, 3) > 2
    error('Error: Input matrices comprise more than two channels')
else
    inchans = size(specTonalLoudness, 3);
    if inchans > 1
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

if ecma == true
    l_start = 1;  % Starting processed signal block for ECMA adherence
else
    l_start = floor(8192/48e3*187.5) + 1;  % Additional term to remove start zero-padding lag
end

%% Define constants

deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2022
c = 0.1618;  % Half-Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2022

dz = 0.5;  % critical band resolution
halfBark = dz:dz:26.5;  % half-critical band rate scale
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2022

% Section 8.1.1 ECMA-418-2:2022
weight_n = 0.5331;  % Equations 113 & 114 ECMA-418-2:2022
% Table 12 ECMA-418-2:2022
a = 0.2918;
b = 0.5459;

% Output sample rate based on tonality hop sizes (Section 6.2.6
% ECMA-418-2:2022)
sampleRate1875 = 48e3/256;

%% Signal processing

% Section 8.1.1 ECMA-418-2:2022
% Weight and combine component specific loudnesses
for chan = inchans:-1:1
    % Equation 114 ECMA-418-2:2022
    maxLoudnessFuncel = a./(max(specTonalLoudness(:, :, chan)...
                                + specNoiseLoudness(:, :, chan), [],...
                                2, "omitnan") + 1e-12) + b;
    specLoudness(:, :, chan) = (specTonalLoudness(:, :, chan).^maxLoudnessFuncel...
                                    + abs((weight_n.*specNoiseLoudness(:, :, chan)).^maxLoudnessFuncel)).^(1./maxLoudnessFuncel);
end

if inchans == 2 && binaural
    % Binaural loudness
    % Section 8.1.5 ECMA-418-2:2022
    specLoudness(:, :, 3) = sqrt(sum(specLoudness.^2, 3)/2);  % Equation 118 
    outchans = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Binaural"];
else
    outchans = inchans;  % assign number of output channels
end

% Section 8.1.2 ECMA-418-2:2022
% Time-averaged specific loudness Equation 115
specLoudnessPowAvg = (sum(specLoudness((59 - l_start):(end + 1 - l_start), :, :).^(1/log10(2)), 1)./size(specLoudness((59 - l_start):(end + 1 - l_start), :, :), 1)).^log10(2);

% Section 8.1.3 ECMA-418-2:2022
% Time-dependent loudness Equation 116
% Discard singleton dimensions
if outchans == 1
    loudnessTDep = sum(specLoudness.*dz, 2);
else
    loudnessTDep = squeeze(sum(specLoudness.*dz, 2));
    specLoudnessPowAvg = squeeze(specLoudnessPowAvg);
end

% Section 8.1.4 ECMA-418-2:2022
% Overall loudness Equation 117
loudnessPowAvg = (sum(loudnessTDep((59 - l_start):(end + 1 - l_start), :).^(1/log10(2)), 1)./size(loudnessTDep((59 - l_start):(end + 1 - l_start), :), 1)).^log10(2);

%% Output assignment

% Assign outputs to structure
if outchans == 3
    loudnessHMS.specLoudness = specLoudness(:, :, 1:2);
    loudnessHMS.specLoudnessPowAvg = specLoudnessPowAvg(:, 1:2);
    loudnessHMS.loudnessTDep = loudnessTDep(:, 1:2);
    loudnessHMS.loudnessPowAvg = loudnessPowAvg(1:2);
    loudnessHMS.specLoudnessBin = specLoudness(:, :, 3);
    loudnessHMS.specLoudnessPowAvgBin = specLoudnessPowAvg(:, 3);
    loudnessHMS.loudnessTDepBin = loudnessTDep(:, 3);
    loudnessHMS.loudnessPowAvgBin = loudnessPowAvg(:, 3);
    loudnessHMS.bandCentreFreqs = bandCentreFreqs;
else
    loudnessHMS.specLoudness = specLoudness;
    loudnessHMS.specLoudnessPowAvg = specLoudnessPowAvg;
    loudnessHMS.loudnessTDep = loudnessTDep;
    loudnessHMS.loudnessPowAvg = loudnessPowAvg;
    loudnessHMS.bandCentreFreqs = bandCentreFreqs;
end

%% Output plotting

% time (s) corresponding with results output
t = (0:(size(specLoudness, 1) - 1))/sampleRate1875;

if outplot
    % Plot figures
    % ------------
    for chan = outchans:-1:1
        cmap_viridis = load('cmap_viridis.txt');
        % Plot results
        fig = figure;
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, t, bandCentreFreqs, permute(specLoudness(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [t(1), t(end) + (t(2) - t(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        %ax1.CLim = [0, ceil(max(specificLoudness(:, :, chan), [], 'all')*10)/10];
        ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
        ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                          "8k", "16k"];
        ax1.YScale = 'log';
        ax1.YLabel.String = 'Frequency, Hz';
        ax1.XLabel.String = 'Time, s';
        ax1.FontName =  'Arial';
        ax1.FontSize = 12;
        colormap(cmap_viridis);
        h = colorbar;
        set(get(h,'label'),'string', {'Specific Loudness,'; 'sone_{HMS}/Bark_{HMS}'});        
        chan_lab = chans(chan);
        title(strcat(chan_lab, ' signal'),...
                     'FontWeight', 'normal', 'FontName', 'Arial');

        ax2 = nexttile(2);
        plot(ax2, t, loudnessPowAvg(1, chan)*ones(size(t)), 'color',...
             cmap_viridis(34, :), 'LineWidth', 0.75, 'DisplayName', "Power" + string(newline) + "time-avg");
        hold on
        plot(ax2, t, loudnessTDep(:, chan), 'color', cmap_viridis(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax2.XLim = [t(1), t(end) + (t(2) - t(1))];
        ax2.YLim = [0, 1.01*ceil(max(loudnessTDep(:, chan))*10)/10];
        ax2.XLabel.String = 'Time, s';
        ax2.YLabel.String = 'Loudness, sone_{HMS}';
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.FontName = 'Arial';
        ax2.FontSize = 12;
        lgd = legend('Location', 'eastoutside', 'FontSize', 8);
        lgd.Title.String = "Overall";
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting if outplot true

% function end
