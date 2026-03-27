function sharpness = acousticSharpFromQuasiLoud(loudQZTDep, specQZLoudness, adjustLoud, fLim, timeStep, sharpWeight, outPlot, binaural)
% sharpness = acousticSharpFromQuasiLoud(loudQZTDep, specQZLoudness,
%                                        adjustLoud, fLim, timeStep,
%                                        sharpWeight, outPlot, binaural)
%
% Returns quasi-sharpness values using quasi-loudness results obtained using
% acousticQuasiLoudZwicker.m or acousticQuasiLoudZwickerWav.m.
%
% The sharpness model used can be specified using the 'sharpMethod' input
% argument. Options comprise 'aures', 'vonBismarck', or 'widmann' (which
% is the model standardised in DIN 45692:2009).
%
% The input argument 'adjustLoud' requires any adjustments applied to the
% quasi-loudness calculations to be declared.
% Optional modifications available comprise adjustments to the spectral
% levels to mimic the application of the ISO 532-3:2023 or
% ECMA-418-2:2025 outer-middle ear filter responses (in 1/3-octaves)
% instead of the ISO 532-1:2017 outer ear transmission, and
% approximated non-linear transformations more closely following
% ISO 532-3 or ECMA-418-2 instead of the ISO 532-1 original.
%
% Since the input matrices will have been calculated using a given sound
% field option ('freeFrontal' or 'diffuse') for the outer ear filter, this
% information is not known to the function, so cannot be included in the
% output.
%
% Inputs
% ------
%
% loudQZTDep : vector or matrix
%   Time-dependent overall quasi-loudness values
%   arranged as [time(, chans)].
%
% specQZLoudness : matrix
%   Specific quasi-loudness values arranged as [time, bands(, chans)].
%
% adjustLoud : keyword string
%   indicates whether adjustments were applied to the loudness for the
%   outer-middle ear filter response and loudness transformations from
%   ISO 532-3:2023 ('iso5323') or ECMA-418-2:2025 ('ecma4182').
%   The option 'iso5321' applies no adjustments, so follows ISO
%   532-1 more closely (for closer agreement with Zwicker's model).
%
% fLim : vector (default: [25, 12600])
%   The frequency limits for the input spectra (these are
%   automatically matched to nearest exact band centre-frequencies).
%
% timeStep : number (default: 0.1)
%   the time step value used for time-dependent inputs.
%
% sharpWeight : keyword string (default: 'aures')
%   the sharpness weighting method to apply. Options: 'aures',
%   'vonbismarck', 'widmann'.
%
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output.
%
% binaural : Boolean true/false (default: true)
%   Flag indicating whether to output binaural sharpness for
%   stereo input signal. (It is assumed the relationship for binaural
%   sharpness follows that of binaural loudness, which seems to be
%   supported by available evidence https://doi.org/10.1051/aacus/2025048)
%
% Returns
% -------
%
% sharpness : structure
%   Contains the output.
%
% sharpness contains the following outputs:
%
% sharpTDep : vector or matrix
%   Time-dependent sharpness arranged as [time(, channels)].
% 
% sharpPowAvg : number or vector
%   Time-power-averaged sharpness arranged as [sharpness(, channels)].
%
% sharp5pcEx : number or vector
%   95th percentile (5% exceeded) sharpness
%   arranged as [sharpness(, channels)].
%
% timeOut : vector
%   Time (seconds) corresponding with time-dependent outputs.
%
% sharpWeight : string
%   Indicates which sharpness weighting method was applied.
%
% adjustLoud : string
%   Indicates which loudness model standard was adjusted for (=adjustLoud
%   input).
%
% If outplot=true, a set of plots is returned illustrating the
% time-dependent sharpness, with the time-aggregated values.
% A set of plots is returned for each input channel.
%
% Assumptions
% -----------
% The input matrices are quasi-Zwicker time-dependent specific and overall
% loudness, with dimensions orientated as
% [Bark bands, time blocks, signal channels] and [time blocks, signal channels].
% The critical band rates for the input specific loudness are the standard
% ISO 532-1 range 0.1 - 24 in 0.1 Bark steps (240 Bark steps).
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
% von Bismarck-Zwicker model is described in:
% 
% von Bismarck, G., 1974. Sharpness as an attribute of the timbre of
% steady sounds, Acta Acustica united with Acustica, 30(3) 159-172.
% https://www.ingentaconnect.com/content/dav/aaua/1974/00000030/00000003/art00006
%
% and:
%
% Zwicker, E. & Fastl, H., 1999. Psychoacoustics: Facts and models.
% Springer-Verlag.
%
% Widmann model is described in DIN 45692:2009
%
% Zwicker loudness is defined in ISO 532-1:2017. Modifications are based on
% ISO 532-3:2023 and ECMA-418-2:2025.
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
% Date created: 30/04/2025
% Date last modified: 27/03/2025
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
        loudQZTDep (:, :) double {mustBeReal}
        specQZLoudness (:, :, :) double {mustBeReal}
        adjustLoud (1, :) string {mustBeMember(adjustLoud,...
                                               {'iso5321',...
                                                'iso5323',...
                                                'ecma4182'})}
        fLim (1, 2) double {mustBeInRange(fLim, 25, 12600)} = [25, 12600]
        timeStep (1, 1) double {mustBePositive} = 0.1
        sharpWeight (1, :) string {mustBeMember(sharpWeight,...
                                                {'aures',...
                                                 'vonbismarck',...
                                                 'widmann'})} = 'aures'
        outPlot (1, 1) {mustBeNumericOrLogical} = false
        binaural (1, 1) {mustBeNumericOrLogical} = true
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks

% Check the inputs are compatible
if size(loudQZTDep, 1) ~= size(specQZLoudness, 1)
    error("Specific and overall loudness inputs must have the same number of time-steps on axis 1.")
end

if ndims(specQZLoudness) - ndims(loudQZTDep) ~= 1 && ...
        ndims(specQZLoudness) - ndims(loudQZTDep) ~= 0
    error("Specific and overall loudness inputs must have compatible dimensions.")
end

% Check the channel number of the input data
if size(specQZLoudness, 3) > 2
    error('Error: Input matrices comprise more than two channels')
else
    chansIn = size(specQZLoudness, 3);
    if chansIn > 1
        chanLabs = ["Stereo left";
                    "Stereo right"];
    else
        chanLabs = "Mono";
    end
end

% fractional octave frequencies
[freqInMid, ~, ~, ~] = noctf(fLim, 3);
[fmAll, ~, ~, ~] = noctf([25, 12600], 3);

%% Define constants

% frequency band limiting indices
fmAllI1 = find(min(freqInMid) == fmAll);
fmAllI2 = find(max(freqInMid) == fmAll);

dz = 0.1;  % Bark number step for ISO 532-1 specific loudness

% ISO 532-1 Table A.8: critical band number (rate) - addition of 1e-4 is as
% as per Annex A.4
barkN = [0.9, 1.8, 2.8, 3.5, 4.4, 5.4, 6.6, 7.9, 9.2, 10.6, 12.3,...
         13.8, 15.2, 16.7, 18.1, 19.3, 20.6, 21.8, 22.7, 23.6, 24] + 1e-4;

% mapped to full frequency band range (omitting top band used only for
% slope)
barkNMap = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1.8, 1.8, 1.8, 2.8, 2.8, 3.5,...
            4.4, 5.4, 6.6, 7.9, 9.2, 10.6, 12.3, 13.8, 15.2, 16.7, 18.1,...
            19.3, 20.6, 21.8, 22.7, 23.6] + 1e-4;

% calculate number of critical bands for output
barkNMapOut = barkNMap(fmAllI1:fmAllI2);
barkNOut = unique(barkNMapOut);
barkI1 = find(min(barkNOut) == barkN);
barkI2 = find(max(barkNOut) == barkN);

if barkI1 == 1
    barkMin = 0.1;
else
    barkMin = barkN(barkI1 - 1) + dz;
end
barkMax = max(round(barkN(barkI2 + 1), 1));
barkAxis = barkMin:dz:barkMax;

% calibration value required to ensure an overall sharpness of 1 acum
% corresponds with a 60 dB (free-field) narrowband noise 1 critical
% bandwidth centred on 1 kHz

switch sharpWeight
    case 'aures'

        switch adjustLoud
            case 'iso5321'
                acum = "A | Q-Z";
                calS = 1.107960224433488;
            case 'iso5323'
                acum = "A | Q-M-G-S";
                calS = 1.098469586031850;
            case 'ecma4182'
                acum = "A | Q-SHM";
                calS = 1.135322210673135;
        end

    case 'vonbismarck'

        switch adjustLoud
            case 'iso5321'
                acum = "vBZ | Q-Z";
                calS = 1.021277699739412;
            case 'iso5323'
                acum = "vBZ | Q-M-G-S";
                calS = 1.021150568220931;
            case 'ecma4182'
                acum = "vBZ | Q-SHM";
                calS = 1.031655615448746;
        end

        q1 = 15;
        q2 = 0.2;
        q3 = 0.308;
        q4 = 0.8;

    case 'widmann'

        switch adjustLoud
            case 'iso5321'
                acum = "W | Q-Z";
                calS = 1.021379909046495;
            case 'iso5323'
                acum = "W | Q-M-G-S";
                calS = 1.021279512048133;
            case 'ecma4182'
                acum = "W | Q-SHM";
                calS = 1.031856782355899;
        end

        q1 = 15.8;
        q2 = 0.15;
        q3 = 0.42;
        q4 = 0.85;
end

%% Signal processing

switch sharpWeight
    case 'aures'
        % Aures weighting
        weightSharp = permute(0.078*(exp(0.171.*permute(repmat(barkAxis, chansIn, 1),...
                                                        [2, 1])))...
                              ./(log(0.05*permute(repmat(loudQZTDep, 1, 1,...
                                                         length(barkAxis)), [3, 2, 1])...
                                     + 1) + eps), [3, 1, 2]);

        % time-dependent sharpness
        % Note: no multiplication by z (otherwise weighting would need /bark term)
        sharpTDep = calS*0.11.*sum(specQZLoudness.*weightSharp*dz, 2);

    case {'vonbismarck', 'widmann'}
        % von Bismarck or Widmann weightings
        weightSharp = ones(1, length(barkAxis));
        weightSharp(barkAxis >= q1) = q2*exp(q3*(barkAxis(barkAxis >= q1) - q1) ) + q4;

        % time-dependent sharpness
        sharpTDep = calS*0.11.*squeeze(sum(specQZLoudness.*weightSharp.*barkAxis*dz, 2))./(loudQZTDep + eps);
end

% Discard singleton dimensions
if chansIn == 2
    sharpTDep = squeeze(sharpTDep);
    chansOut = chansIn;

    % Binaural sharpness
    if binaural
        chansOut = 3;
        sharpTDepBin = sqrt(sum(sharpTDep(:, 1:2).^2, 2)/2);
        sharpPowAvgBin = (sum(sharpTDepBin.^(1/log10(2)), 1)./size(sharpTDepBin, 1)).^log10(2);
        sharp5pcExBin = prctile(sharpTDepBin, 95, 1);
        chanLabs = [chanLabs; "Binaural"];
    end
else
    chansOut = chansIn;
end

% ensure any 0 loudness values also have 0 sharpness
sharpTDep(loudQZTDep == 0) = 0;

% overall (power-averaged) sharpness
sharpPowAvg = (sum(sharpTDep.^(1/log10(2)), 1)./size(sharpTDep, 1)).^log10(2);

% overall (95th percentile) sharpness
sharp5pcEx = prctile(sharpTDep, 95, 1);

% time (s) corresponding with results output
timeOut = (0:(size(specQZLoudness, 1) - 1))*timeStep;

%% Output assignment

% Assign outputs to structure
sharpness.sharpTDep = sharpTDep;
sharpness.sharpPowAvg = sharpPowAvg;
sharpness.sharp5pcEx = sharp5pcEx;
sharpness.timeOut = timeOut.';
sharpness.sharpWeight = sharpWeight;
sharpness.adjustLoud = adjustLoud;

if binaural && chansIn == 2
    sharpness.sharpTDepBin = sharpTDepBin;
    sharpness.sharpPowAvgBin = sharpPowAvgBin;
    sharpness.sharp5pcExBin = sharp5pcExBin;
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
        plot(ax, timeOut, sharpPowAvg(1, chan)*ones(size(timeOut)), 'color',...
             cmap_viridis(34, :), 'LineWidth', 1, 'DisplayName', "Power" + string(newline) + "time-avg");
        hold on
        plot(ax, timeOut, sharp5pcEx(1, chan)*ones(size(timeOut)), 'color',...
             cmap_viridis(34, :), 'LineWidth', 1, 'LineStyle', ':',...
             'DisplayName', "5%" + string(newline) + "exceeded");
        plot(ax, timeOut, sharpTDep(:, chan), 'color', cmap_viridis(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        if max(sharpTDep(:, chan)) > 0
            ax.YLim = [0, 1.1*ceil(max(sharpTDep(:, chan))*10)/10];
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
        chan_lab = chanLabs(chan);
        title(strcat(chan_lab, ' signal'),...
                     'FontWeight', 'normal', 'FontName', 'Arial');
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting

end % end of acousticSharpFromQuasiLoud function
