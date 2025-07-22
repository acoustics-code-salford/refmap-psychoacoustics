function loudnessSHM = acousticSHMLoudness(p, sampleRateIn, axisN, soundField, waitBar, outPlot, binaural)
% loudnessSHM = acousticSHMLoudness(p, sampleRateIn, axisN, soundField,
%                                   outPlot, binaural)
%
% Returns loudness values according to ECMA-418-2:2025 (using the Sottek
% Hearing Model) for an input calibrated single mono or single stereo
% audio (sound pressure) time-series signal, p. For stereo signals, the
% binaural loudness can be calculated, and each channel is also analysed
% separately.
%
% Inputs
% ------
% p : vector or 2D matrix
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% sampleRateIn : integer
%                the sample rate (frequency) of the input signal(s)
%
% axisN : integer (1 or 2, default: 1)
%         the time axis along which to calculate the tonality
%
% waitBar : keyword string (default: true)
%           determines whether a progress bar displays during processing
%           (set waitBar to false for doing multi-file parallel calculations)
%
% soundField : keyword string (default: 'freeFrontal')
%              determines whether the 'freeFrontal' or 'diffuse' field stages
%              are applied in the outer-middle ear filter, or 'noOuter' uses
%              only the middle ear stage, or 'noEar' omits ear filtering.
%              Note: these last two options are beyond the scope of the
%              standard, but may be useful if recordings made using
%              artificial outer/middle ear are to be processed using the
%              specific recorded responses.
%
% outPlot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% binaural : Boolean true/false (default: true)
%            flag indicating whether to output combined binaural loudness
%            for stereo input signal.
% 
% Returns
% -------
%
% loudnessSHM : structure
%               contains the output
%
% loudnessSHM contains the following outputs:
%
% specLoudness : matrix
%                time-dependent specific loudness for each critical band
%                arranged as [time, bands(, channels)]
%
% specloudnessPowAvg : matrix
%                      time-power-averaged specific loudness for each
%                      critical band
%                      arranged as [bands(, channels)]
%
% specTonalLoudness : matrix
%                     time-dependent specific tonal loudness for each
%                     (half) critical band
%                     arranged as [time, bands(, channels)]
%
% specNoiseLoudness : matrix
%                     time-dependent specific noise loudness for each
%                     (half) critical band
%                     arranged as [time, bands(, channels)]
%
% loudnessTDep : vector or matrix
%                 time-dependent overall loudness
%                 arranged as [time(, channels)]
% 
% loudnessPowAvg : number or vector
%                  time-power-averaged overall loudness
%                  arranged as [loudness(, channels)]
%
% bandCentreFreqs : vector
%                   centre frequencies corresponding with each critical
%                   band rate
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% soundField : string
%              identifies the soundfield type applied (= input argument)
%
% If binaural=true, a corresponding set of outputs for the combined
% binaural loudness are also contained in loudnessSHM
%
% If outplot=true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall loudness, with the latter also indicating the time-aggregated
% value. A set of plots is returned for each input channel, with another
% set for the combined binaural loudness, if binaural=true. In that case,
% the indicated sound level corresponds with the channel with the highest
% sound level.
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in Pascals
% (Pa).
%
% Requirements
% ------------
% Signal Processing Toolbox
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 22/09/2023
% Date last modified: 22/07/2025
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
% This code calls sub-component file 'cmap_viridis.txt'. The contents of
% the file includes a copy of data obtained from the repository 
% https://github.com/BIDS/colormap, and is CC0 1.0 licensed for modified
% use, see https://creativecommons.org/publicdomain/zero/1.0 for
% information.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        p (:, :) double {mustBeReal}
        sampleRateIn (1, 1) double {mustBePositive, mustBeInteger}
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        soundField (1, :) string {mustBeMember(soundField,...
                                               {'freeFrontal',...
                                                'diffuse', ...
                                                'noOuter', ...
                                                'noEar'})} = 'freeFrontal'
        waitBar {mustBeNumericOrLogical} = true
        outPlot {mustBeNumericOrLogical} = false
        binaural {mustBeNumericOrLogical} = true
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks
% Orient input matrix
if axisN == 2
    p = p.';
end

% Check the length of the input data (must be longer than 300 ms)
if size(p, 1) <  300/1000*sampleRateIn
    error('Error: Input signal is too short along the specified axis to calculate loudness (must be longer than 300 ms)')
end

% Check the channel number of the input data
if size(p, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    chansIn = size(p, 2);
    if chansIn > 1
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

%% Define constants

sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2025 [r_s]
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2025 [deltaf(f=0)]
c = 0.1618;  % Half-overlapping Bark band centre-frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2025

dz = 0.5;  % critical band overlap [deltaz]
halfBark = 0.5:dz:26.5;  % half-overlapping critical band rate scale [z]
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025 [F(z)]

% Section 8.1.1 ECMA-418-2:2025
weight_n = 0.5331;  % Equations 113 & 114 ECMA-418-2:2025 [w_n]
% Table 12 ECMA-418-2:2025
a = 0.2918;
b = 0.5459;

% Output sample rate based on tonality hop sizes (Section 6.2.6
% ECMA-418-2:2025) [r_sd]
sampleRate1875 = sampleRate48k/256;

% standardised epsilon
epsilon = 1e-12;

%% Signal processing

% Input pre-processing
% --------------------
if sampleRateIn ~= sampleRate48k  % Resample signal
    [p_re, ~] = shmResample(p, sampleRateIn);
else  % don't resample
    p_re = p;
end

% Calculate specific loudnesses for tonal and noise components
% ------------------------------------------------------------

% Obtain tonal and noise component specific loudnesses from Sections 5 & 6 ECMA-418-2:2025
tonalitySHM = acousticSHMTonality(p_re, sampleRate48k, 1, soundField, waitBar,...
                                  false);

specTonalLoudness = tonalitySHM.specTonalLoudness;  % [N'_tonal(l,z)]
specNoiseLoudness = tonalitySHM.specNoiseLoudness;  % [N'_noise(l,z)]

% Section 8.1.1 ECMA-418-2:2025
% Weight and combine component specific loudnesses
for chan = chansIn:-1:1

    % Equation 114 ECMA-418-2:2025 [e(z)]
    maxLoudnessFuncel = a./(max(specTonalLoudness(:, :, chan)...
                                + specNoiseLoudness(:, :, chan), [],...
                                2, "omitnan") + epsilon) + b;

    % Equation 113 ECMA-418-2:2025 [N'(l,z)]
    specLoudness(:, :, chan) = (specTonalLoudness(:, :, chan).^maxLoudnessFuncel...
                                + abs((weight_n.*specNoiseLoudness(:, :, chan)).^maxLoudnessFuncel)).^(1./maxLoudnessFuncel);
end

if chansIn == 2 && binaural
    % Binaural loudness
    % Section 8.1.5 ECMA-418-2:2025 Equation 118 [N'_B(l,z)]
    specLoudness(:, :, 3) = sqrt(sum(specLoudness(:, :, 1:2).^2, 3)/2);
    specTonalLoudness(:, :, 3) = sqrt(sum(specTonalLoudness(:, :, 1:2).^2, 3)/2);
    specNoiseLoudness(:, :, 3) = sqrt(sum(specNoiseLoudness(:, :, 1:2).^2, 3)/2);
    chansOut = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Binaural"];
else
    chansOut = chansIn;  % assign number of output channels
end

% Section 8.1.2 ECMA-418-2:2025
% Time-averaged specific loudness Equation 115 [N'(z)]
specLoudnessPowAvg = (sum(specLoudness((57 + 1):end, :, :).^(1/log10(2)), 1)./size(specLoudness((57 + 1):end, :, :), 1)).^log10(2);

% Section 8.1.3 ECMA-418-2:2025
% Time-dependent loudness Equation 116 [N(l)]
% Discard singleton dimensions
if chansOut == 1
    loudnessTDep = sum(specLoudness.*dz, 2);
    specLoudnessPowAvg = transpose(specLoudnessPowAvg);
else
    loudnessTDep = squeeze(sum(specLoudness.*dz, 2));
    specLoudnessPowAvg = squeeze(specLoudnessPowAvg);
end

% Section 8.1.4 ECMA-418-2:2025
% Overall loudness Equation 117 [N]
loudnessPowAvg = (sum(loudnessTDep((57 + 1):end, :).^(1/log10(2)), 1)./size(loudnessTDep((57 + 1):end, :), 1)).^log10(2);

% time (s) corresponding with results output [t]
timeOut = (0:(size(specLoudness, 1) - 1))/sampleRate1875;

%% Output plotting

if outPlot
    % Plot figures
    % ------------
    for chan = chansOut:-1:1
        cmap_viridis = load('cmap_viridis.txt');
        % Plot results
        fig = figure;
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, timeOut, bandCentreFreqs, permute(specLoudness(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        ax1.CLim = [0, ceil(max(specLoudness(:, :, chan), [], 'all')*10)/10];
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
        set(get(h,'label'),'string', {'Specific Loudness,'; 'sone_{SHM}/Bark_{SHM}'});        
        chan_lab = chans(chan);

        % Create A-weighting filter
        weightFilt = weightingFilter('A-weighting', sampleRate48k);
        % Filter signal to determine A-weighted time-averaged level
        if chan == 3
            pA = weightFilt(p_re);
            LAeq2 = 20*log10(rms(pA, 1)/2e-5);
            % take the higher channel level as representative (PD ISO/TS
            % 12913-3:2019 Annex D)
            [LAeq, LR] = max(LAeq2);
            % if branch to identify which channel is higher
            if LR == 1
                whichEar = ' left ear';
            else
                whichEar = ' right ear';
            end  % end of if branch

            chan_lab = chan_lab + whichEar;

        else
            pA = weightFilt(p_re(:, chan));
            LAeq = 20*log10(rms(pA)/2e-5);
        end
        
        title(strcat(chan_lab,...
                     ' signal sound pressure level =', {' '},...
                     num2str(round(LAeq,1)), "dB {\itL}_{Aeq}"),...
                     'FontWeight', 'normal', 'FontName', 'Arial');

        ax2 = nexttile(2);
        plot(ax2, timeOut, loudnessPowAvg(1, chan)*ones(size(timeOut)), 'color',...
             cmap_viridis(34, :), 'LineWidth', 1, 'DisplayName', "Power" + string(newline) + "time-avg");
        hold on
        plot(ax2, timeOut, loudnessTDep(:, chan), 'color', cmap_viridis(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        if max(loudnessTDep(:, chan)) > 0
            ax2.YLim = [0, 1.1*ceil(max(loudnessTDep(:, chan))*10)/10];
        end
        ax2.XLabel.String = 'Time, s';
        ax2.YLabel.String = 'Loudness, sone_{SHM}';
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.GridAlpha = 0.075;
        ax2.GridLineStyle = '--';
        ax2.GridLineWidth = 0.25;
        ax2.FontName = 'Arial';
        ax2.FontSize = 12;
        lgd = legend('Location', 'eastoutside', 'FontSize', 8);
        lgd.Title.String = "Overall";
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting if outplot true

%% Output assignment

% Assign outputs to structure
if chansOut == 3
    loudnessSHM.specLoudness = specLoudness(:, :, 1:2);
    loudnessSHM.specTonalLoudness = specTonalLoudness(:, :, 1:2);
    loudnessSHM.specNoiseLoudness = specNoiseLoudness(:, :, 1:2);
    loudnessSHM.specLoudnessPowAvg = specLoudnessPowAvg(:, 1:2);
    loudnessSHM.loudnessTDep = loudnessTDep(:, 1:2);
    loudnessSHM.loudnessPowAvg = loudnessPowAvg(1:2);
    loudnessSHM.specLoudnessBin = specLoudness(:, :, 3);
    loudnessSHM.specTonalLoudnessBin = specTonalLoudness(:, :, 3);
    loudnessSHM.specNoiseLoudnessBin = specNoiseLoudness(:, :, 3);
    loudnessSHM.specLoudnessPowAvgBin = specLoudnessPowAvg(:, 3);
    loudnessSHM.loudnessTDepBin = loudnessTDep(:, 3);
    loudnessSHM.loudnessPowAvgBin = loudnessPowAvg(:, 3);
    loudnessSHM.bandCentreFreqs = bandCentreFreqs;
    loudnessSHM.timeOut = timeOut;
    loudnessSHM.soundField = soundField;
else
    loudnessSHM.specLoudness = specLoudness;
    loudnessSHM.specTonalLoudness = specTonalLoudness;
    loudnessSHM.specNoiseLoudness = specNoiseLoudness;
    loudnessSHM.specLoudnessPowAvg = specLoudnessPowAvg;
    loudnessSHM.loudnessTDep = loudnessTDep;
    loudnessSHM.loudnessPowAvg = loudnessPowAvg;
    loudnessSHM.bandCentreFreqs = bandCentreFreqs;
    loudnessSHM.timeOut = timeOut;
    loudnessSHM.soundField = soundField;
end

% end of function
