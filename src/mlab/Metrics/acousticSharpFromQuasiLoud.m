function sharpness = acousticSharpFromQuasiLoud(loudQZTDep, specQZLoudness, timeStep, method, outPlot)
% sharpness = acousticSharpFromQuasiLoud(loudQZTDep, specQZLoudness,
%                                             method, outPlot)
%
% Returns quasi-sharpness values using quasi-loudness results obtained using
% acousticQuasiLoudZwicker.m or acousticQuasiLoudZwickerWav.m.
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
%
% loudQZTDep : vector or matrix
%              the time-dependent overall loudness values
%              arranged as [time(, chans)]
%
% specQZLoudness : matrix
%                  the specific loudness values
%                  arranged as [time, bands(, chans)]
%
% timeStep : number
%            the time step value used for time-dependent inputs
%
% method : keyword string (default: 'aures')
%          the sharpness method to apply. Options: 'aures', 'vonbismarck',
%          'widmann'
%
% outPlot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% Returns
% -------
%
% sharpness : structure
%             contains the output
%
% sharpness contains the following outputs:
%
% sharpTDep : vector or matrix
%                 time-dependent sharpness
%                 arranged as [time(, channels)]
% 
% sharpPowAvg : number or vector
%                   time-power-averaged sharpness
%                   arranged as [sharpness(, channels)]
%
% sharp5pcEx : number or vector
%                  95th percentile (5% exceeded) sharpness
%                  arranged as [sharpness(, channels)]
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% method : string
%          indicates which sharpness method was applied
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
% von Bismarck model is described in:
% 
% von Bismarck, G., 1974. Sharpness as an attribute of the timbre of
% steady sounds, Acta Acustica united with Acustica, 30(3) 159-172.
% https://www.ingentaconnect.com/content/dav/aaua/1974/00000030/00000003/art00006
%
% Widmann model is described in DIN 45692:2009
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
        loudQZTDep (:, :) double {mustBeReal}
        specQZLoudness (:, :, :) double {mustBeReal}
        timeStep (1, 1) double {mustBePositive}
        method (1, :) string {mustBeMember(method,...
                                           {'aures',...
                                            'vonbismarck',...
                                            'widmann'})} = 'aures'
        outPlot {mustBeNumericOrLogical} = false
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
    chans = size(specQZLoudness, 3);
    if chans > 1
        chanLabs = ["Stereo left";
                 "Stereo right"];
    else
        chanLabs = "Mono";
    end
end

%% Define constants

dz = 0.1;  % Bark number step for ISO 532-1 specific loudness
bark = dz:dz:24;  % Bark numbers for ISO 532-1 specific loudness

% calibration value required to ensure an overall sharpness of 1 acum
% corresponds with a 60 dB (free-field) narrowband noise 1 critical
% bandwidth centred on 1 kHz

switch method
    case 'aures'
        calS = 0.971827112506741;
        acum = "Aur | Quasi-Zwicker";

    case 'vonbismarck'
        calS = 0.972308112404852;
        acum = "vBis | Quasi-Zwicker";
        q1 = 15;
        q2 = 0.2;
        q3 = 0.308;
        q4 = 0.8;

    case 'widmann'
        calS = 0.975032094599035;
        acum = "Widm | Quasi-Zwicker";
        q1 = 15.8;
        q2 = 0.15;
        q3 = 0.42;
        q4 = 0.85;

end

%% Signal processing

switch method
    case 'aures'
        % Aures weighting
        weightSharp = permute(0.078*(exp(0.171.*permute(repmat(bark, chans, 1),...
                                                        [2, 1])))...
                              ./(log(0.05*permute(repmat(loudQZTDep, 1, 1,...
                                                         length(bark)), [3, 2, 1])...
                                     + 1) + eps), [3, 1, 2]);

        % time-dependent sharpness
        % Note: no multiplication by z (otherwise weighting would need /bark term)
        sharpTDep = calS*0.11.*sum(specQZLoudness.*weightSharp*dz, 2);

    case {'vonbismarck', 'widmann'}
        % von Bismarck or Widmann weightings
        weightSharp = ones(1, length(bark));
        weightSharp(bark>=q1) = q2*exp(q3*(bark(bark>=q1) - q1) ) + q4;

        % time-dependent sharpness
        sharpTDep = calS*0.11.*squeeze(sum(specQZLoudness.*weightSharp.*bark*dz, 2))./(loudQZTDep + eps);

end

% Discard singleton dimensions
if chans ~= 1
    sharpTDep = squeeze(sharpTDep);
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
sharpness.method = method;

%% Output plotting

if outPlot
    % Plot figures
    % ------------
    for chan = chans:-1:1
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
end  % end of if branch for plotting if outplot true

end % end of acousticSharpFromQuasiLoud function
