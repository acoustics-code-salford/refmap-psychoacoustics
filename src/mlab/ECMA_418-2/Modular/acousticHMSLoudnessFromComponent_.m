function [loudnessPowAvg, loudnessTimeVar, specificLoudness,...
          specificLoudnessPowAvg, bandCentreFreqs]...
          = acousticHMSLoudnessFromComponent_(specificTonalLoudness,...
                                              specificNoiseLoudness,...
                                              outplot, ecma, binaural)
% [loudnessPowAvg, loudnessTimeVar, specificLoudness, 
%  specificLoudnessPowAvg, bandCentreFreqs]
%  = acousticHMSLoudnessFromComponent(specificTonalLoudness,
%                                     specificNoiseLoudness,
%                                     outplot, ecma, binaural)
%
% Returns loudness values according to ECMA-418-2:2022 (using the Hearing
% Model of Sottek) for input component specific tonal loudness and specific
% noise loudness, obtained using acousticHMSTonality().
%
% Inputs
% ------
% specificTonalLoudness : 2D or 3D matrix
%                         the specific tonal loudness values calculated for
%                         a sound pressure signal (single mono or single
%                         stereo audio)
%
% specificNoiseLoudness : 2D or 3D matrix
%                         the specific noise loudness values calculated for
%                         a sound pressure signal (single mono or single
%                         stereo audio)
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% ecma : Boolean true/false (default: true)
%        flag indicating whether to maintain strict standard adherence to
%        ECMA-418-2:2022 Equation 40, or otherwise to use an alternative
%        that provides closer time-alignment of the time-dependent tonality
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
% Date last modified: 15/03/2024
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
        specificTonalLoudness (:, :, :) double {mustBeReal}
        specificNoiseLoudness (:, :, :) double {mustBeReal}
        outplot {mustBeNumericOrLogical} = false
        ecma {mustBeNumericOrLogical} = true
        binaural {mustBeNumericOrLogical} = true
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))

%% Input checks
% Check the size of the input matrices (must match)
if ~isequal(size(specificTonalLoudness), size(specificNoiseLoudness))
    error('Error: Input loudness matrix sizes must match')
end

% Check the channel number of the input data
if size(specificTonalLoudness, 3) > 2
    error('Error: Input matrices comprise more than two channels')
else
    inchans = size(specificTonalLoudness, 3);
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
c = 0.1618;  % Half-Bark band centre-frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2022

halfBark = 0.5:0.5:26.5;  % half-critical band rate scale
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2022

% Section 8.1.1 ECMA-418-2:2022
weight_n = 0.5331;  % Equations 113 & 114 ECMA-418-2:2022
% Table 12 ECMA-418-2:2022
a = 0.2918;
b = 0.5459;


%% Signal processing

% Section 8.1.1 ECMA-418-2:2022
% Weight and combine component specific loudnesses
for chan = inchans:-1:1
    % Equation 114 ECMA-418-2:2022
    maxLoudnessFuncel = a./(max(specificTonalLoudness(:, :, chan)...
                                + specificNoiseLoudness(:, :, chan), [],...
                                1, "omitnan") + 1e-12) + b;
    specificLoudness(:, :, chan) = (specificTonalLoudness(:, :, chan).^maxLoudnessFuncel...
                                    + abs((weight_n.*specificNoiseLoudness(:, :, chan)).^maxLoudnessFuncel)).^(1./maxLoudnessFuncel);
end

if inchans > 1 && binaural
    % Binaural loudness
    % Section 8.1.5 ECMA-418-2:2022
    specificLoudness = sqrt(sum(specificLoudness.^2, 3)/2);  % Equation 118 
    outchans = 1;  % output number of 'channels' to single binaural
else
    outchans = inchans;  % assign number of output channels
end

% Section 8.1.2 ECMA-418-2:2022
% Time-averaged specific loudness Equation 115
specificLoudnessPowAvg = (sum(specificLoudness(:, (59 - l_start):(end + 1 - l_start), :).^(1/log10(2)), 2)./size(specificLoudness(:, (59 - l_start):(end + 1 - l_start), :), 2)).^log10(2);

% Section 8.1.3 ECMA-418-2:2022
% Time-dependent loudness Equation 116
if outchans == 1
    loudnessTimeVar = shiftdim(sum(specificLoudness.*0.5, 1), 1);
else
    loudnessTimeVar = squeeze(sum(specificLoudness.*0.5, 1));
end

% Section 8.1.4 ECMA-418-2:2022
% Overall loudness Equation 117
loudnessPowAvg = (sum(loudnessTimeVar((59 - l_start):(end + 1 - l_start), :).^(1/log10(2)), 1)./size(loudnessTimeVar((59 - l_start):(end + 1 - l_start), :), 1)).^log10(2);

%% Output plotting

% time (s) corresponding with results output
t = (0:(size(specificLoudness, 2) - 1))/187.5;

if outplot
    % Plot figures
    % ------------
    for chan = outchans:-1:1
        cmap_viridis = load('cmap_viridis.txt');
        % Plot results
        if inchans > 1 && binaural
            chan_lab = "Binaural";
        else
            chan_lab = chans(chan);
        end
        fig = figure;
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, t, bandCentreFreqs, specificLoudness(:, :, chan),...
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
        
        title(strcat(chan_lab, ' signal'),...
                     'FontWeight', 'normal', 'FontName', 'Arial');

        ax2 = nexttile(2);
        plot(ax2, t, loudnessPowAvg(1, chan)*ones(size(t)), 'color',...
             [0.1, 0.3, 0.9], 'LineWidth', 0.75, 'DisplayName', "Power" + string(newline) + "time-avg");
        hold on
        plot(ax2, t, loudnessTimeVar(:, chan), 'color', [0.1, 0.9, 0.6],...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "varying");
        hold off
        ax2.XLim = [t(1), t(end) + (t(2) - t(1))];
        ax2.YLim = [0, 1.01*ceil(max(loudnessTimeVar(:, chan))*10)/10];
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
