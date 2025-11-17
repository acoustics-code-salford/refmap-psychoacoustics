function annoyance = acousticPsychAnnoy(p, sampleRateIn, axisN, startSkip, soundField, outPlot)
% annoyance = acousticPsychAnnoy(p, sampleRateIn, axisN, startSkip, soundField, prctExceed, outPlot)
%
% Returns psychoacoustic annoyance values according to Widmann (1992) [also
% known as 'Zwicker Psychoacoustic Annoyance - see below].
% for an input calibrated single mono or single stereo audio
% (sound pressure) time-series signal, p. For stereo signals, the binaural
% roughness can be calculated, and each channel is also analysed separately.
%
% Widmann, U. (1992). Ein Modell der Psychoakustischen Lästigkeit von
% Schallen und seine Anwendung in der Praxis der Lärmbeurteilung (A model
% of the psychoacoustic annoyance of sounds and its application in noise
% assessment practice) [Doctoral thesis, Technische Universität München
% (Technical University of Munich)].
%
% Commonly misattributed to Zwicker, E. and Fastl, H. Second ed.,
% Psychoacoustics, Facts and Models, 2nd ed. Springer-Verlag, Berlin, 1999.
%
% Lotinga, M. J. B. and A. J. Torija (2025). "Comment on "A study on
% calibration methods of noise annoyance data from listening tests"
% [J. Acoust. Soc. Am. 156, 1877–1886 (2024)]." Journal of the Acoustical
% Society of America 157(5): 3282–3285.
%
% Inputs
% ------
% p : vector or 2D matrix
%   Input signal as single mono or stereo audio (sound pressure) signals.
%
% sampleRateIn : integer
%   Sample rate (frequency) of the input signal(s).
%
% axisN : integer (1 or 2, default: 1)
%   Time axis along which to calculate the tonality.
%
% startSkip : number (default: 2)
%   Amount of time to skip at the start of the signal for
%   calculating time-aggregated outputs (starts from next input
%   sample). The default value is obtained from the result of
%   calculating fluctuation strength for the relevant reference
%   signal.
%
% soundField : keyword string (default: 'freeFrontal')
%   Determines whether the 'freeFrontal' or 'diffuse' field stages
%   are applied in the outer-middle ear filtering.
%
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output.
%
% Returns
% -------
%
% annoyance : structure
%   contains the output
%
% annoyance contains the following outputs:
%
% loudness : structure
%   contains sub-outputs:
%
%   NTDep : vector or 2D matrix
%       Time-dependent overall loudness arranged as [time(, channels)].
%
%   specN : matrix
%       Time-dependent specific loudness for each critical band
%       arranged as [time, bands(, channels)].
%
%   N5 : number or vector
%       Time-aggregated loudness, as 5%-exceeded (95th percentile).
%
%   NPwAvg : number or vector
%       Time-aggregated loudness, as power-average (see ECMA-418-2).
%
%   timeOutN : vector
%       Time (seconds) corresponding with time-dependent loudness outputs.
%
% sharpness : structure
%   contains sub-outputs:
%
%   STDep : vector or 2D matrix
%       Time-dependent sharpness arranged as [time(, channels)].
%
%   S5 : number or vector
%       Time-aggregated sharpness, as 5%-exceeded (95th percentile).
%
%   SPwAvg : number or vector
%       Time-aggregated sharpness, as power-average (see ECMA-418-2).
%
%   timeOutS : vector
%       Time (seconds) corresponding with time-dependent sharpness outputs.
%
% fluctuation : structure
%   contains sub-outputs:
%
%   FTDep : vector or 2D matrix
%       Time-dependent overall fluctuation strength arranged as
%       [time(, channels)].
%
%   specF : matrix
%       Time-dependent specific fluctuation strength for each critical band
%       arranged as [time, bands(, channels)].
%
%   F5 : number or vector
%       Time-aggregated fluctuation strength, as 5%-exceeded (95th percentile).
%
%   F10 : number or vector
%       Time-aggregated fluctuation strength, as 10%-exceeded (90th
%       percentile).
%
%   timeOutF : vector
%       Time (seconds) corresponding with time-dependent fluctuation 
%       strength outputs.
%
% roughness : structure
%   contains sub-outputs:
%
%   roughnessTDep : vector or 2D matrix
%       Time-dependent overall roughness arranged as [time(, channels)].
%
%   specRoughness : matrix
%       Time-dependent specific roughness for each critical band
%       arranged as [time, bands(, channels)].
%
%   R5 : number or vector
%       Time-aggregated roughness, as 5%-exceeded (95th percentile).
%
%   timeOutR : vector
%       Time (seconds) corresponding with time-dependent roughness outputs.
%
% psychAnnoyTDep : vector or 2D matrix
%   Time-dependent psychoacoustic annoyance arranged as [time(, channels)].
%
% psychAnnoyFrom5Pc : number or vector
%   Overall psychoacoustic annoyance based on 5%-exceeded time-aggregated
%   metric values.
%
% psychAnnoyFromAlt : number or vector
%   Overall psychoacoustic annoyance based on alternative time-aggregated
%   metric values (thought to be better estimates of overall perception for
%   time-varying sound qualities).
%
% timeOutPA : vector
%   Time (seconds) corresponding with time-dependent psychacoustic annoyance output.
%
% soundField : string
%              identifies the soundfield type applied (= input argument)
%
% If outplot=true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall roughness, with the latter also indicating the time-aggregated
% value. A set of plots is returned for each input channel.
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
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 07/07/2025
% Date last modified: 17/11/2025
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
% This code calls sub-component file 'cmap_inferno.txt'. The contents of
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
        startSkip (1, 1) {mustBePositive} = 2     
        soundField (1, :) string {mustBeMember(soundField,...
                                               {'freeFrontal',...
                                                'diffuse'})} = 'freeFrontal'
        outPlot (1, 1) {mustBeNumericOrLogical} = false
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))


%% Input checks
% Orient input matrix
if axisN == 2
    p = p.';
end

% Check the length of the input data (must be longer than 2.5 s)
if size(p, 1) <=  2.5*sampleRateIn
    error("Error: Input signal is too short along the specified axis to calculate psychoacoustic annoyance (must be longer than 2.5 s due to fluctuation strength calculation, which skips at least 2 s from start and 0.5 s from end of output calculation)")
end

% Check the channel number of the input data
if size(p, 2) > 2
    error("Error: Input signal comprises more than two channels")
else
    chansIn = size(p, 2);
    if chansIn == 2
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

% convert soundField keyword for compatibility with Audio Toolbox
if strcmp(soundField, 'freeFrontal')
    soundFieldKey = "free";
else
    soundFieldKey = "diffuse";
end

%% Define constants

% calibration value to ensure reference signal equals 1 sone (reference:
% sinusoid at 1 kHz at 40 dB free-field)
calN = 0.999511042933208;
% calibration value to ensure reference signal equals 1 asper (reference:
% sinusoid at 1 kHz at 60 dB free-field, 100% amplitude modulated at 70 Hz
% modulation)
calR = 1.073321839408406;
% calibration value to ensure reference signal equals 1 vacil (reference:
% sinusoid at 1 kHz at 60 dB free-field, 100% amplitude modulated at 4 Hz
% modulation)
calF = 1.044096149945446;
% calibration value to ensure reference signal equals 1 acum (reference:
% narrowband noise 1 critical bandwidth centred at 1 kHz at 60 dB
% free-field - Zwicker's reference rather than DIN 45692)
calS = 0.993681909452517;

% sampling periods for output component metrics
samplePerN = 1/500;  % period for loudness
samplePerS = 1/500;  % period for sharpness
samplePerF = 1/500;  % period for fluctuation strength
samplePerR = 1/2000;  % period for roughness

barkAxis = 0.5:0.5:23.5;  % critical band rates

% sample rates
sampleRateN = 1/samplePerN;  % rate for loudness
sampleRateS = 1/samplePerS;  % rate for loudness
sampleRateF = 1/samplePerF;  % rate for loudness
sampleRateR = 1/samplePerR;  % rate for roughness
sampleRateMax = max([sampleRateN, sampleRateS, sampleRateF, sampleRateR]);
samplePerMax = 1/sampleRateMax;

% end skip in seconds, to avoid processing artefacts for modulation metrics
endSkip = 0.5;

%% Signal processing

[loudness, specLoudness] = acousticLoudness(p, sampleRateIn, 1,...
                                            'SoundField', soundFieldKey,...
                                            'Method', 'ISO 532-1',...
                                            'TimeVarying', true);
loudness = calN*loudness;
specLoudness = calN*specLoudness;

% calculation of other sound qualities (note: although specific loudness
% can be used to calculate R and F, the resolutions required differ, so
% that loudness would need to be calculated from pressure twice,
% diminishing the efficiency benefit)
[roughness, specRoughness] = acousticRoughness(p, sampleRateIn, 1,...
                                               'SoundField', soundFieldKey);
roughness = calR*roughness;
specRoughness = calR*specRoughness;

[fluctuation, specFluctuation] = acousticFluctuation(specLoudness);
fluctuation = calF*fluctuation;
specFluctuation = calF*specFluctuation;

sharpness = acousticSharpness(specLoudness, 'TimeVarying', true,...
                                   'Weighting', 'DIN 45692');
sharpness = calS*sharpness;

% calculate weighted sharpness component
weightSharp = zeros(size(sharpness));
weightSharp(sharpness > 1.75) = 0.25*log10(loudness(sharpness > 1.75) + 10)...
                       .*(sharpness(sharpness > 1.75) - 1.75);

% interpolation to match metric sampling rates
lenN = size(loudness, 1);
lenS = size(weightSharp, 1);
lenF = size(fluctuation, 1);
lenR = size(roughness, 1);
lenMax = max([lenN, lenS, lenF, lenR]);
timeOutN = 0:samplePerN:lenN*samplePerN - samplePerN;
timeOutS = 0:samplePerS:lenS*samplePerS - samplePerS;
timeOutF = 0:samplePerF:lenF*samplePerF - samplePerF;
timeOutR = 0:samplePerR:lenR*samplePerR - samplePerR;
timeOutMax = 0:samplePerMax:lenMax*samplePerMax - samplePerMax;

if lenN < lenMax
    loudUp = interp1(timeOutN, loudness, timeOutMax, 'pchip').';
end

if lenS < lenMax
    weightSharpUp = interp1(timeOutS, weightSharp, timeOutMax, 'pchip').';
end

if lenF < lenMax
    fluctUp = interp1(timeOutF, fluctuation, timeOutMax, 'pchip').';
end

if lenR < lenMax
    roughUp = interp1(timeOutR, roughness, timeOutMax, 'pchip').';
end

% weighted time-dependent modulation
weightMod = 2.18./loudUp.^0.4.*(0.4*fluctUp + 0.6*roughUp);

% time-dependent psychoacoustic annoyance
psychAnnoyTDep = loudUp.*(1 + sqrt(weightSharpUp.^2 + weightMod.^2));

% convert skips to output samples for time aggregation
startSkipNSF = startSkip*sampleRateNSF + 1;
endSkipNSF = endSkip*sampleRateNSF;
startSkipR = startSkip*sampleRateR + 1;
endSkipR = endSkip*sampleRateR;

% time-aggregated psychoacoustic annoyance from 5%-exceeded values
N5 = prctile(loudness(startSkipNSF:end - endSkipNSF, :), 95, 1);
S5 = prctile(sharpness(startSkipNSF:end - endSkipNSF, :), 95, 1);
F5 = prctile(fluctuation(startSkipNSF:end - endSkipNSF, :), 95, 1);
R5 = prctile(roughness(startSkipR:end - endSkipR, :), 95, 1);
psychAnnoyFrom5Pc = PA(N5, S5, F5, R5);

% time-aggregated psychoacoustic annoyance from alternative values
NPwAvg = mean(loudness(startSkipNSF:end - endSkipNSF, :).^(1/log10(2)), 1).^(log10(2));
SPwAvg = mean(sharpness(startSkipNSF:end - endSkipNSF, :).^(1/log10(2)), 1).^(log10(2));
F10 = prctile(fluctuation(startSkipNSF:end - endSkipNSF, :), 90, 1);
R10 = prctile(roughness(startSkipR:end - endSkipR, :), 90, 1);
psychAnnoyFromAlt = PA(NPwAvg, SPwAvg, F10, R10);

timeOutPA = 0:samplePerMax:size(psychAnnoyTDep, 1)*samplePerMax - samplePerMax;

%% Output plotting

if outPlot
    % TODO needs completing (based on detection)
    % cmap_plasma = load('cmap_plasma.txt');
    % cmap_inferno = load('cmap_inferno.txt');
    % cmap_magma = load('cmap_magma.txt');
    % cmap_viridis = load('cmap_viridis.txt');
    % cmap_cividis = load('cmap_cividis.txt');
    % cmap_mako = load('cmap_mako.txt');
    % 
    % % Plot figures
    % % ------------
    % for chan = chansIn:-1:1
    %     % Plot results
    %     fig = figure;
    %     set(fig, 'position', [300, 200, 1300, 600])
    %     tiledlayout(fig, 2, 3);
    %     movegui(fig, 'center');
    %     ax1 = nexttile(1);
    % 
    % 
    %     fig = figure;
    %     set(fig, 'position', [300, 200, 1300, 600])
    %     tiledlayout(fig, 2, 3);
    %     movegui(fig, 'center');
    % 
    %     ax1 = nexttile(1);
    %     surf(ax1, [timeOutN, timeOutN(end) + timeOutN(2)] - timeOutN(2)/2, fm_vis - fBandWidth_vis/2,...
    %          arrangeSpectro(dBSpecTarget(:, :, targChan)),...
    %          EdgeColor='none', FaceColor='interp');
    %     view(2);
    %     set(ax1, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
    %         'YScale', 'log', 'YLim', [f1(1), f2(end)],...
    %         'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
    %         'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
    %     ax1.YLabel.String = "Frequency, Hz";
    %     ax1.XLabel.String = "Time, s";
    %     colormap(ax1, cmap_viridis);
    %     cbar = colorbar; clim([0, max(dBSpecTarget(:, :, targChan), [], 'all')]);
    %     cbar.Label.String = "dB re 2e-5 Pa";
    %     ax1.Title.String = "Target spectrogram";
    %     ax1.TitleFontWeight = "normal";
    % 
    %     ax2 = nexttile(2);
    %     surf(ax2, [t, t(end) + timeStep] - timeStep/2, fm_vis - fBandWidth_vis/2,...
    %          arrangeSpectro(detectabilitydB(:, :, targChan)),...
    %          EdgeColor='none', FaceColor='interp');
    %     view(2);
    %     set(ax2, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
    %         'YScale', 'log', 'YLim', [f1(1), f2(end)],...
    %         'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
    %         'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
    %     ax2.YLabel.String = "Frequency, Hz";
    %     ax2.XLabel.String = "Time, s";
    %     colormap(ax2, cmap_magma);
    %     cbar = colorbar; clim([0, max(detectabilitydB(:, :, targChan), [], 'all')]);
    %     cbar.Label.String = "Detectability 10log_{10}{\it d'}, dB";
    %     ax2.Title.String = "Masked target detectability";
    %     ax2.TitleFontWeight = "normal";
    % 
    %     ax3 = nexttile(3);
    %     plot(ax3, t, dBATDepTarget(:, targChan), 'color', cmap_magma(34, :), 'DisplayName', "Target")
    %     hold on
    %     plot(ax3, t, dBATDepMasker(:, maskChan), 'color', cmap_magma(166, :), 'DisplayName', "Masker")
    %     plot(ax3, t, dBATDepDiscTarget(:, targChan), 'color', cmap_magma(100, :),...
    %         'LineStyle', ':', 'LineWidth', 2, 'DisplayName', "Target discounted")
    %     hold off
    %     set(ax3, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
    %         'YLim', [min(dBATDepMasker(:, maskChan), [], 'all') - 10,...
    %                  max(max(dBATDepTarget(:, targChan), [], 'all'),...
    %                  max(dBATDepMasker(:, maskChan), [], 'all'))*1.05],...
    %         'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.075,...
    %         'GridLineStyle', '--', 'GridLineWidth', 0.25);
    %     ax3.YLabel.String = "dB(A) re 2e-5 Pa";
    %     ax3.XLabel.String = "Time, s";
    %     ax3.Title.String = "Time-dependent levels";
    %     ax3.TitleFontWeight = "normal";
    %     legend(ax3, 'Location', 'eastoutside');
    % 
    %     ax4 = nexttile(4);
    %     surf(ax4, [t, t(end) + timeStep] - timeStep/2, fm_vis - fBandWidth_vis/2,...
    %          arrangeSpectro(dBSpecMasker(:, :, targChan)),...
    %          EdgeColor='none', FaceColor='interp');
    %     view(2);
    %     set(ax4, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
    %         'YScale', 'log', 'YLim', [f1(1), f2(end)],...
    %         'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
    %         'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
    %     ax4.YLabel.String = "Frequency, Hz";
    %     ax4.XLabel.String = "Time, s";
    %     colormap(ax4, cmap_viridis);
    %     cbar = colorbar; clim([0, max(dBSpecMasker(:, :, targChan), [], 'all')]);
    %     cbar.Label.String = "dB re 2e-5 Pa";
    %     ax4.Title.String = "Masker spectrogram";
    %     ax4.TitleFontWeight = "normal";
    % 
    %     ax5 = nexttile(5);
    %     surf(ax5, [t, t(end) + timeStep] - timeStep/2, fm_vis - fBandWidth_vis/2,...
    %          arrangeSpectro(dBSpecDiscTarget(:, :, targChan)),...
    %          EdgeColor='none', FaceColor='interp');
    %     view(2);
    %     set(gca, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
    %         'YScale', 'log', 'YLim', [f1(1), f2(end)],...
    %         'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
    %         'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
    %     ax5.YLabel.String = "Frequency, Hz";
    %     ax5.XLabel.String = "Time, s";
    %     colormap(ax5, cmap_magma);
    %     cbar = colorbar; clim([0, max(dBSpecDiscTarget(:, :, targChan), [], 'all')]);
    %     cbar.Label.String = "Detectability-discounted \newline     level, dB re 2e-5 Pa";
    %     ax5.Title.String = "Target detectability-discounted spectrogram";
    %     ax5.TitleFontWeight = "normal";
    % 
    %     ax6 = nexttile(6);
    %     levelVals = [LAETarget(targChan), LAEMasker(maskChan), LAEDiscTarget(targChan)];
    %     labelVals = num2cell(round(levelVals, 1));
    %     labelCats = ["Target", "Masker", "Discount. targ."];
    %     % a trick using stacked bar to get the legend mapped
    %     b = bar(ax6,...
    %             diag(levelVals, 0),...
    %             'stacked', 'FaceColor', 'flat');
    %     set(b, {'FaceColor'}, {cmap_magma(34, :); cmap_magma(166, :); cmap_magma(100, :)})
    %     set(ax6, 'YLim', [min(levelVals, [], 'all')/1.25, max(levelVals)*1.1])
    %     ax6.YLabel.String = "{\it L}_{AE}, dB re 2e-5 Pa";
    %     ax6.XTickLabel = [];
    %     lg6 = legend(ax6, labelCats, 'Location','eastoutside');
    %     lg6.Direction = 'normal';
    %     % data labels
    %     y = sum(reshape(cell2mat(get(b', 'YData')),size(b, 2), []), 1); 
    %     x = unique(cell2mat(get(b', 'XData')),'stable');
    %     offset = range(ylim)*.1; 
    %     text(x, y - offset, labelVals, 'HorizontalAlignment', 'Center',...
    %          'VerticalAlignment', 'bottom', 'Color', 'w');
    %     ax6.Title.String = "Overall levels";
    %     ax6.TitleFontWeight = "normal";

    % end
end

%% Output assignment

annoyance.loudness.NTDep = loudness;
annoyance.loudness.specN = specLoudness;
annoyance.loudness.N5 = N5;
annoyance.loudness.NPwAvg = NPwAvg;
annoyance.loudness.timeOutN = timeOutN;

annoyance.sharpnessSTDep = sharpness;
annoyance.sharpness.S5 = S5;
annoyance.sharpness.SPwAvg = SPwAvg;
annoyance.loudness.timeOutS = timeOutS;

annoyance.fluctuation.FTDep = fluctuation;
annoyance.fluctuation.specF = specFluctuation;
annoyance.fluctuation.F5 = F5;
annoyance.fluctuation.F10 = F10;
annoyance.loudness.timeOutF = timeOutF;

annoyance.roughness.RTDep = roughness;
annoyance.roughness.specR = specRoughness;
annoyance.roughness.R5 = R5;
annoyance.roughness.R10 = R10;
annoyance.loudness.timeOutR = timeOutR;

annoyance.psychAnnoyTDep = psychAnnoyTDep;
annoyance.psychAnnoyFrom5Pc = psychAnnoyFrom5Pc;
annoyance.psychAnnoyFromAlt = psychAnnoyFromAlt;
annoyance.timeOutPA = timeOutPA;

annoyance.soundField = soundField;

%% nested function to calculate psychoacoustic annoyance
function pa = PA(N, S, F, R)
    wS = zeros(size(S));
    wS(S > 1.75) = 0.25*log10(N(S > 1.75) + 10).*(S(S > 1.75) - 1.75);
    wFR = 2.18./N.^0.4.*(0.4*F + 0.6*R);

    pa = N.*(1 + sqrt(wS.^2 + wFR.^2));
end  % end of nested function for psychoacoustic annoyance

end

