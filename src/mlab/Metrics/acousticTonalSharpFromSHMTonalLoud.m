function tonalSharpnessSHM = acousticTonalSharpFromSHMTonalLoud(specSHMTonalLoudness, method, outPlot, binaural)
% tonalSharpnessSHM = acousticTonalSharpFromSHMTonalLoud(specSHMTonalLoudness, method,
%                                                        outPlot, binaural)
%
% Returns tonal sharpness values using Sottek Hearing Model (SHM) specific
% tonal loudness obtained using acousticSHMLoudness.m or 
% acousticSHMTonality.m. This is faster than calculating via
% acousticTonalSharpnessSHM.m (which calls acousticSHMLoudness.m).
%
% The sharpness model used can be specified using the 'method' input
% argument. Options comprise 'aures', 'vonbismarck', or 'widmann' (which
% is the model standardised in DIN 45692:2009).
%
% Since the input matrices will have been calculated using a given sound
% field option for the outer ear filter, this information is not known to
% the function, so cannot be included in the output.
%
% Inputs
% ------
% specSHMTonalLoudness : matrix
%   The specific SHM tonal loudness values calculated for
%   a sound pressure signal (single mono or single
%   stereo audio) arranged as [time, bands(, chans)]
%
% method : keyword string (default: 'aures')
%   The sharpness method to apply. Options: 'aures', 'vonbismarck',
%   'widmann'
%
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output
%
% binaural : Boolean true/false (default: true)
%   Flag indicating whether to output binaural tonal sharpness for
%   stereo input signal. (It is assumed the relationship for binaural tonal
%   sharpness follows that of binaural SHM loudness, which seems to be
%   supported by available evidence https://doi.org/10.1051/aacus/2025048)
% 
% Returns
% -------
%
% tonalSharpnessSHM : structure
%   Contains the output
%
% tonalSharpnessSHM contains the following outputs:
%
% tonalSharpnessTDep : vector or matrix
%   Time-dependent tonal sharpness arranged as [time(, channels)]
% 
% tonalSharpnessPowAvg : number or vector
%   Time-power-averaged tonal sharpness
%   arranged as [tonalSharpness(, channels)]
%
% timeOut : vector
%   Time (seconds) corresponding with time-dependent outputs
%
% method : string
%   Indicates which sharpness method was applied
%
% If binaural=true, a corresponding set of outputs for the binaural
% tonal sharpness are also contained in tonalSharpnessSHM
%
% If outPlot=true, a plot is returned illustrating the
% time-dependent tonal sharpness, with the time-aggregated value.
% A plot is returned for each input channel, with another
% plot for the binaural tonal sharpness, if binaural=true.
%
% Assumptions
% -----------
% The input matrices are ECMA-418-2:2025 specific tonal loudness, with
% dimensions orientated as [bands, time blocks, signal channels].
% The critical band rates for the input specific loudness are the standard
% ECMA-418-2 range 0.5 - 26.5 in 0.5 Bark steps (53 half-bands).
% The time step for the input specific loudness is 1/187.5.
%
% References
% ----------
% Aures model is described in:
%
% Aures, W., 1985. Berechnungsverfahren für den sensorischen Wohlklang
% beliebiger Schallsignale, Acta Acustica united with Acustica, 59(2),
% 130-141.
% https://www.ingentaconnect.com/content/dav/aaua/1985/00000059/00000002/art00008
%
% von Bismarck model is described in:
% 
% von Bismarck, G., 1974. Sharpness as an attribute of the timbre of
% steady sounds, Acta Acustica united with Acustica, 30(3) 159-172.
% https://www.ingentaconnect.com/content/dav/aaua/1974/00000030/00000003/art00006
%
% Widmann model is described in DIN 45692:2009
%
% The implementation here incorporates the basic approach to modifying the
% model for different loudness inputs according to:
%
% Swift, S.H. & Gee, K.L., 2017. Extending sharpness calculation for an alternative
% loudness metric input, Journal of the Acoustical Society of America, 142(6), EL549–EL554.
% https://doi.org/10.1121/1.5016193
%
% The modification is undertaken using an expression for deriving the Bark
% critical band rate from frequency according to:
%
% Volk, F., 2015. Comparison and fitting of analytical expressions to
% existing data for the critical-band concept, Acta Acustica united with
% Acustica, 101(6), 1157-1167. https://doi.org/10.3813/AAA.918908
%
% The assumed binaural perception of sharpness is based on evidence found
% in:
%
% Hochbaum, F, Hundt,  T, Fiebig, A & Brinkmann, F, 2025. Directional
% sharpness perception under different listening conditions, Acta Acustica,
% 9, 60. https://doi.org/10.1051/aacus/2025048
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
% Date created: 11/12/2024
% Date last modified: 11/12/2025
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
        specSHMTonalLoudness (:, :, :) double {mustBeReal}
        method (1, :) string {mustBeMember(method,...
                                           {'aures',...
                                            'vonbismarck',...
                                            'widmann'})} = 'aures'
        outPlot {mustBeNumericOrLogical} = false
        binaural {mustBeNumericOrLogical} = true
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks

% Check the channel number of the input data
if size(specSHMTonalLoudness, 3) > 2
    error('Error: Input matrices comprise more than two channels')
else
    chansIn = size(specSHMTonalLoudness, 3);
    if chansIn > 1
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

%% Define constants

deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2025
c = 0.1618;  % Half-Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2025
dz = 0.5;  % critical band resolution [deltaz]
halfBark = dz:dz:26.5;  % half-critical band rate scale
f = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025
% z = 13*atan(0.76*(f/1000)) + 3.5*atan((f/7500).^2);  % Bark Eq 6.1 Fastl
% & Zwicker (no longer used)
z = 32.12*(1 - (1 + (f/873.47).^1.18).^-0.4);  % Bark eq 9 Volk, 2015

dt = 1/187.5;  % time step (resolution, s)

% calibration value required to ensure an overall sharpness of 1 acum
% corresponds with a 60 dB (free-field) sinusoid at 1 kHz

switch method
    case 'aures'
        calS = 1.050805334544524;
        acum = "Aur | SHM";

    case 'vonbismarck'
        calS = 0.983541704055661;
        acum = "vBis | SHM";
        q1 = 15;
        q2 = 0.2;
        q3 = 0.308;
        q4 = 0.8;

    case 'widmann'
        calS = 0.984912403132019;
        acum = "Widm | SHM";
        q1 = 15.8;
        q2 = 0.15;
        q3 = 0.42;
        q4 = 0.85;

end

%% Signal processing

if chansIn == 2 && binaural
    % Binaural loudness
    % Section 8.1.5 ECMA-418-2:2025 Equation 118
    specSHMTonalLoudness(:, :, 3) = sqrt(sum(specSHMTonalLoudness.^2, 3)/2);
    chansOut = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Binaural"];
else
    chansOut = chansIn;  % assign number of output channels
end

% Time-dependent loudness
% Discard singleton dimensions
if chansOut == 1
    tonalLoudnessTDep = sum(specSHMTonalLoudness.*dz, 2);
else
    tonalLoudnessTDep = squeeze(sum(specSHMTonalLoudness.*dz, 2));
end

switch method
    case 'aures'
        % Aures weighting
        weightSharp = permute(0.078*(exp(0.171.*permute(repmat(z, chansOut, 1),...
                                                        [2, 1])))...
                              ./(log(0.05*permute(repmat(tonalLoudnessTDep, 1, 1,...
                                                         length(z)), [3, 2, 1])...
                                     + 1) + eps), [3, 1, 2]);
        % Adjustment to the first part (1-57 time steps) to avoid visualisation of
        % nonsense values generated from first redundant 0.3s of time-dependent
        % loudness
        weightSharp(1:57, :, :) = repmat(weightSharp(58, :, :), 57, 1, 1);

        % time-dependent sharpness
        % Note: no multiplication by z (otherwise weighting would need /z term)
        tonalSharpnessTDep = calS*0.11.*sum(specSHMTonalLoudness.*weightSharp*dz, 2);

    case {'vonbismarck', 'widmann'}
        % von Bismarck or Widmann weightings
        weightSharp = ones(1, length(z));
        weightSharp(z>=q1) = q2*exp(q3*(z(z>=q1) - q1) ) + q4;

        % time-dependent sharpness
        tonalSharpnessTDep = calS*0.11.*squeeze(sum(specSHMTonalLoudness.*weightSharp.*z*dz, 2))./(tonalLoudnessTDep + eps);
        % Adjustment to the first part (1-57 time steps) to avoid visualisation of
        % nonsense values generated from first redundant 0.3s of time-dependent
        % loudness
        tonalSharpnessTDep(1:57, :, :) = calS*0.11.*squeeze(sum(specSHMTonalLoudness(1:57, :, :).*weightSharp.*z*dz, 2))./tonalLoudnessTDep(58, :, :);

end

% Discard singleton dimensions
if chansOut ~= 1
    tonalSharpnessTDep = squeeze(tonalSharpnessTDep);
end

% ensure any 0 loudness values also have 0 sharpness
tonalSharpnessTDep(tonalLoudnessTDep==0) = 0;

% overall (power-averaged) sharpness
tonalSharpnessPowAvg = (sum(tonalSharpnessTDep((57 + 1):end, :).^(1/log10(2)), 1)./size(tonalSharpnessTDep((57 + 1):end, :), 1)).^log10(2);

if chansIn == 2 && binaural
    % alternative binaural calculation method (Hochbaum et al, 2025)
    tonalSharpnessTDepBin2 = sqrt(sum(tonalSharpnessTDep.^2, 2)/2);
    tonalSharpnessPowAvgBin2 = (sum(tonalSharpnessTDepBin2((57 + 1):end, :).^(1/log10(2)), 1)./size(tonalSharpnessTDepBin2((57 + 1):end, :), 1)).^log10(2);
end

% time (s) corresponding with results output
timeOut = transpose((0:(size(specSHMTonalLoudness, 1) - 1))*dt);

%% Output assignment

% Assign outputs to structure
if chansOut == 3
    tonalSharpnessSHM.tonalSharpnessTDep = tonalSharpnessTDep(:, 1:2);
    tonalSharpnessSHM.tonalSharpnessPowAvg = tonalSharpnessPowAvg(1:2);
    tonalSharpnessSHM.tonalSharpnessTDepBin = tonalSharpnessTDep(:, 3);
    tonalSharpnessSHM.tonalSharpnessPowAvgBin = tonalSharpnessPowAvg(:, 3);
    tonalSharpnessSHM.tonalSharpnessTDepBin2 = tonalSharpnessTDepBin2;
    tonalSharpnessSHM.tonalSharpnessPowAvgBin2 = tonalSharpnessPowAvgBin2;
    tonalSharpnessSHM.timeOut = timeOut;
    tonalSharpnessSHM.method = method;
else
    tonalSharpnessSHM.tonalSharpnessTDep = tonalSharpnessTDep;
    tonalSharpnessSHM.tonalSharpnessPowAvg = tonalSharpnessPowAvg;
    tonalSharpnessSHM.timeOut = timeOut;
    tonalSharpnessSHM.method = method;
end

%% Output plotting

if outPlot
    % Plot figures
    % ------------
    for chan = chansOut:-1:1
        cmap_viridis = load('cmap_viridis.txt');
        % Plot results
        fig = figure;
        movegui(fig, 'center');
        ax = gca();
        plot(ax, timeOut, tonalSharpnessPowAvg(1, chan)*ones(size(timeOut)), ':', 'color',...
             cmap_viridis(34, :), 'LineWidth', 1.5, 'DisplayName', "Power" + string(newline) + "time-avg");
        hold on
        plot(ax, timeOut, tonalSharpnessTDep(:, chan), 'color', cmap_viridis(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        if max(tonalSharpnessTDep(:, chan)) > 0
            ax.YLim = [0, 1.1*ceil(max(tonalSharpnessTDep(:, chan))*10)/10];
        end
        ax.XLabel.String = "Time, s";
        ax.YLabel.String = "Tonal sharpness, acum_{{\it T}, " + acum + "}";
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridAlpha = 0.075;
        ax.GridLineStyle = '--';
        ax.GridLineWidth = 0.25;
        ax.FontName = 'Arial';
        ax.FontSize = 12;
        lgd = legend('Location', 'eastoutside', 'FontSize', 8);
        lgd.Title.String = "Overall";
        chan_lab = chans(chan);
        title(strcat(chan_lab, ' signal'),...
                     'FontWeight', 'normal', 'FontName', 'Arial');
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting

% function end
