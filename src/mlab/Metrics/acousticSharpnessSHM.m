function sharpnessSHM = acousticSharpnessSHM(p, sampleRatein, axisn, fieldtype, method, outplot, binaural)
% sharpnessSHM = acousticSharpnessSHM(p, sampleRatein, axisn,...
%                                     fieldtype, method, outplot, binaural)
%
% Returns sharpness values using Sottek Hearing Model (SHM) loudness for an
% input calibrated single mono or single stereo audio (sound pressure)
% time-series signal, p. For stereo signals, the binaural sharpness can
% also be calculated, and each channel also analysed separately.
%
% The sharpness model used can be specified using the 'method' input
% argument. Options comprise 'aures', 'vonBismarck', or 'widmann' (which
% is the model standardised in DIN 45692:2009).
%
% Since the input matrices will have been calculated using a given sound
% field option ('free-frontal' or 'diffuse') for the outer ear filter, this
% information is not known to the function, so cannot be included in the
% output.
%
% Inputs
% ------
% p : vector or 2D matrix
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% sampleRatein : integer
%                the sample rate (frequency) of the input signal(s)
%
% axisn : integer (1 or 2, default: 1)
%         the time axis along which to calculate the sharpness
%
% fieldtype : keyword string (default: 'free-frontal')
%             determines whether the 'free-frontal' or 'diffuse' field stages
%             are applied in the outer-middle ear filter
%
% method : keyword string (default: 'aures')
%          the sharpness method to apply. Options: 'aures', 'vonbismarck',
%          'widmann'
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% binaural : Boolean true/false (default: true)
%            flag indicating whether to output binaural sharpness for
%            stereo input signal. (Experimental: it is assumed the
%            relationship for binaural sharpness follows that of binaural
%            SHM loudness.
% 
% Returns
% -------
%
% sharpnessSHM : structure
%                     contains the output
%
% sharpnessSHM contains the following outputs:
%
% sharpnessTDep : vector or matrix
%                 time-dependent sharpness
%                 arranged as [time(, channels)]
% 
% sharpnessPowAvg : number or vector
%                   time-power-averaged sharpness
%                   arranged as [sharpness(, channels)]
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% soundField : string
%              identifies the soundfield type applied (the input argument
%              fieldtype)
%
% method : string
%          indicates which sharpness method was applied
%
% If binaural=true, a corresponding set of outputs for the binaural
% sharpness are also contained in sharpnessSHM
%
% If outplot=true, a plot is returned illustrating the
% time-dependent sharpness, with the time-aggregated value.
% A plot is returned for each input channel, with another
% plot for the binaural sharpness, if binaural=true.
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in Pascals
% (Pa).
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
% Date created: 01/11/2024
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
        sampleRatein (1, 1) double {mustBePositive, mustBeInteger}
        axisn (1, 1) {mustBeInteger, mustBeInRange(axisn, 1, 2)} = 1
        fieldtype (1, :) string {mustBeMember(fieldtype,...
                                              {'free-frontal',...
                                               'diffuse',...
                                               'noOuter'})} = 'free-frontal'
        method (1, :) string {mustBeMember(method,...
                                           {'aures',...
                                            'vonbismarck',...
                                            'widmann'})} = 'aures'
        outplot {mustBeNumericOrLogical} = false
        binaural {mustBeNumericOrLogical} = true
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks
% Orient input matrix
if axisn == 2
    p = p.';
end

% Check the length of the input data (must be longer than 300 ms)
if size(p, 1) <  300/1000*sampleRatein
    error('Error: Input signal is too short along the specified axis to calculate sharpness (must be longer than 300 ms)')
end

% Check the channel number of the input data
if size(p, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    inchans = size(p, 2);
    if inchans > 1
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

%% Define constants

deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2024
c = 0.1618;  % Half-Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2024
dz = 0.5;  % critical band resolution [deltaz]
halfBark = dz:dz:26.5;  % half-critical band rate scale
f = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2024
% z = 13*atan(0.76*(f/1000)) + 3.5*atan((f/7500).^2);  % Bark Eq 6.1 Fastl
% & Zwicker (no longer used)
z = 32.12*(1 - (1 + (f/873.47).^1.18).^-0.4);  % Bark eq 9 Volk, 2015

dt = 1/187.5;  % time step (resolution, s)

% calibration value required to ensure an overall sharpness of 1 acum
% corresponds with a 60 dB (free-field) narrowband noise 1 critical
% bandwidth centred on 1 kHz

switch method
    case 'aures'
        calS = 1.0643285;
        acum = "Aur | SHM";

    case 'vonbismarck'
        calS = 0.98393;
        acum = "vBis | SHM";
        q1 = 15;
        q2 = 0.2;
        q3 = 0.308;
        q4 = 0.8;

    case 'widmann'
        calS = 0.9854;
        acum = "Widm | SHM";
        q1 = 15.8;
        q2 = 0.15;
        q3 = 0.42;
        q4 = 0.85;

end

%% Signal processing

% Calculate specific loudness
% ---------------------------

% Obtain specific loudness from ECMA-418-2:2024
loudnessSHM = acousticSHMLoudness(p, sampleRatein, 1, fieldtype, false);

specSHMLoudness = loudnessSHM.specLoudness;

if inchans == 2 && binaural
    % Binaural loudness
    % Section 8.1.5 ECMA-418-2:2024 Equation 118
    specSHMLoudness(:, :, 3) = loudnessSHM.specLoudnessBin;
    outchans = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Binaural"];
else
    outchans = inchans;  % assign number of output channels
end

% Time-dependent loudness
% Discard singleton dimensions
if outchans == 1
    loudnessTDep = sum(specSHMLoudness.*dz, 2);
else
    loudnessTDep = squeeze(sum(specSHMLoudness.*dz, 2));
end

switch method
    case 'aures'
        % Aures weighting
        weightSharp = permute(0.078*(exp(0.171.*permute(repmat(z, outchans, 1),...
                                                        [2, 1])))...
                              ./(log(0.05*permute(repmat(loudnessTDep, 1, 1,...
                                                         length(z)), [3, 2, 1])...
                                     + 1) + eps), [3, 1, 2]);
        % Adjustment to the first part (1-57 time steps) to avoid visualisation of
        % nonsense values generated from first redundant 0.3s of time-dependent
        % loudness
        weightSharp(1:57, :, :) = repmat(weightSharp(58, :, :), 57, 1, 1);

        % time-dependent sharpness
        % Note: no multiplication by z (otherwise weighting would need /z term)
        sharpnessTDep = calS*0.11.*sum(specSHMLoudness.*weightSharp*dz, 2);

    case {'vonbismarck', 'widmann'}
        % von Bismarck or Widmann weightings
        weightSharp = ones(1, length(z));
        weightSharp(z>=q1) = q2*exp(q3*(z(z>=q1) - q1) ) + q4;

        % time-dependent sharpness
        sharpnessTDep = calS*0.11.*squeeze(sum(specSHMLoudness.*weightSharp.*z*dz, 2))./(loudnessTDep + eps);
        % Adjustment to the first part (1-57 time steps) to avoid visualisation of
        % nonsense values generated from first redundant 0.3s of time-dependent
        % loudness
        sharpnessTDep(1:57, :, :) = calS*0.11.*squeeze(sum(specSHMLoudness(1:57, :, :).*weightSharp.*z*dz, 2))./loudnessTDep(58, :, :);
end

% Discard singleton dimensions
if outchans ~= 1
    sharpnessTDep = squeeze(sharpnessTDep);
end

% ensure any 0 loudness values also have 0 sharpness
sharpnessTDep(loudnessTDep==0) = 0;

% overall (power-averaged) sharpness
sharpnessPowAvg = (sum(sharpnessTDep((57 + 1):end, :).^(1/log10(2)), 1)./size(sharpnessTDep((57 + 1):end, :), 1)).^log10(2);

% time (s) corresponding with results output
timeOut = (0:(size(specSHMLoudness, 1) - 1))*dt;

%% Output assignment

% Assign outputs to structure
if outchans == 3
    sharpnessSHM.sharpnessTDep = sharpnessTDep(:, 1:2);
    sharpnessSHM.sharpnessPowAvg = sharpnessPowAvg(1:2);
    sharpnessSHM.sharpnessTDepBin = sharpnessTDep(:, 3);
    sharpnessSHM.sharpnessPowAvgBin = sharpnessPowAvg(:, 3);
    sharpnessSHM.timeOut = timeOut;
    sharpnessSHM.soundField = fieldtype;
    sharpnessSHM.method = method;
else
    sharpnessSHM.sharpnessTDep = sharpnessTDep;
    sharpnessSHM.sharpnessPowAvg = sharpnessPowAvg;
    sharpnessSHM.timeOut = timeOut;
    sharpnessSHM.soundField = fieldtype;
    sharpnessSHM.method = method;
end

%% Output plotting

if outplot
    % Plot figures
    % ------------
    for chan = outchans:-1:1
        cmap_viridis = load('cmap_viridis.txt');
        % Plot results
        fig = figure;
        movegui(fig, 'center');
        ax = gca();
        plot(ax, timeOut, sharpnessPowAvg(1, chan)*ones(size(timeOut)), 'color',...
             cmap_viridis(34, :), 'LineWidth', 1, 'DisplayName', "Power" + string(newline) + "time-avg");
        hold on
        plot(ax, timeOut, sharpnessTDep(:, chan), 'color', cmap_viridis(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        if max(sharpnessTDep(:, chan)) > 0
            ax.YLim = [0, 1.1*ceil(max(sharpnessTDep(:, chan))*10)/10];
        end
        ax.XLabel.String = "Time, s";
        ax.YLabel.String = "Sharpness, acum_{" + acum + "}";
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
end  % end of if branch for plotting if outplot true

% function end
