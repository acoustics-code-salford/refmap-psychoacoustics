function annoyance = acousticPsychAnnoy(p, sampleRateIn, axisN, startSkip, soundField, prctExceed, outPlot)
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
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% sampleRateIn : integer
%                the sample rate (frequency) of the input signal(s)
%
% axisN : integer (1 or 2, default: 1)
%         the time axis along which to calculate the tonality
%
% startSkip : number (default: 2)
%             the amount of time to skip at the start of the signal for
%             calculating time-aggregated outputs (starts from next input
%             sample). The default value is obtained from the result of
%             calculating fluctuation strength for the relevant reference
%             signal
%
% soundField : keyword string (default: 'freeFrontal')
%              determines whether the 'freeFrontal' or 'diffuse' field stages
%              are applied in the outer-middle ear filtering.
%
% prctExceed : vector (default: [0, 1, 5, 10, 50, 90, 95, 99, 100]
%              time-aggregation percent-exceeded for calculation
%              (calculated from 100-statistical percentile)
%
% outPlot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% Returns
% -------
%
% annoyance : structure
%                contains the output
%
% annoyance contains the following outputs:
%
% specRoughness : matrix
%                 time-dependent specific roughness for each (half)
%                 critical band
%                 arranged as [time, bands(, channels)]
%
% specRoughnessAvg : matrix
%                    time-averaged specific roughness for each (half)
%                    critical band
%                    arranged as [bands(, channels)]
%
% roughnessTDep : vector or matrix
%                 time-dependent overall roughness
%                 arranged as [time(, channels)]
% 
% roughness90Pc : number or vector
%                 time-aggregated (90th percentile) overall roughness
%                 arranged as [roughness(, channels)]
%
% bandCentreFreqs : vector
%                   centre frequencies corresponding with each (half)
%                   critical band rate scale width
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% soundField : string
%              identifies the soundfield type applied (the input argument
%              fieldtype)
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
% Date last modified: 10/07/2025
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
        prctExceed (1, :) {mustBeReal, mustBeInteger} = [0, 1, 5, 10, 50, 90, 95, 99, 100]
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
samplePerNSF = 1/500;  % period for loudness, sharpness and fluctuation strength
samplePerR = 1/2000;  % period for roughness
sampleRateNSF = 500;  % rate for loudness, sharpness and fluctuation strength
sampleRateR = 2000;  % rate for roughness

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
timeOutNSF = 0:samplePerNSF:size(loudness, 1)*samplePerNSF - samplePerNSF;
timeOutRPA = 0:samplePerR:size(roughness, 1)*samplePerR - samplePerR;
loudUp = interp1(timeOutNSF, loudness, timeOutRPA, 'pchip').';
fluctUp = interp1(timeOutNSF, fluctuation, timeOutRPA, 'pchip').';
weightSharpUp = interp1(timeOutNSF, weightSharp, timeOutRPA, 'pchip').';

weightMod = 2.18./loudUp.^0.4.*(0.4*fluctUp + 0.6*roughness);

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

%% Output plotting

if outPlot
    cmap_plasma = load('cmap_plasma.txt');
    cmap_inferno = load('cmap_inferno.txt');
    cmap_magma = load('cmap_magma.txt');
    cmap_viridis = load('cmap_viridis.txt');
    cmap_cividis = load('cmap_cividis.txt');

    % Plot figures
    % ------------
    for chan = chansIn:-1:1
        % Plot results
        fig = figure;
        set(fig, 'position', [300, 200, 1300, 600])
        tiledlayout(fig, 2, 3);
        movegui(fig, 'center');
        ax1 = nexttile(1);
    end
end

%% Output assignment

annoyance.loudness.loudnessTDep = loudness;
annoyance.loudness.specLoudness = specLoudness;
annoyance.loudness.N5 = N5;
annoyance.loudness.NPwAvg = NPwAvg;

annoyance.sharpness.sharpnessTDep = sharpness;
annoyance.sharpness.S5 = S5;
annoyance.sharpness.SPwAvg = SPwAvg;

annoyance.fluctuation.fluctuationTDep = fluctuation;
annoyance.fluctuation.specFluctuation = specFluctuation;
annoyance.fluctuation.F5 = F5;
annoyance.fluctuation.F10 = F10;

annoyance.roughness.roughnessTDep = roughness;
annoyance.roughness.specRoughness = specRoughness;
annoyance.roughness.R5 = R5;
annoyance.roughness.R10 = R10;

annoyance.psychAnnoyTDep = psychAnnoyTDep;
annoyance.psychAnnoyFrom5Pc = psychAnnoyFrom5Pc;
annoyance.psychAnnoyFromAlt = psychAnnoyFromAlt;

annoyance.timeOut.timeOutNSF = timeOutNSF;
annoyance.timeOut.timeOutRPA = timeOutRPA;

%% nested function to calculate psychoacoustic annoyance
function pa = PA(N, S, F, R)
    wS = zeros(size(S));
    wS(S > 1.75) = 0.25*log10(N(S > 1.75) + 10).*(S(S > 1.75) - 1.75);
    wFR = 2.18./N.^0.4.*(0.4*F + 0.6*R);

    pa = N.*(1 + sqrt(wS.^2 + wFR.^2));
end  % end of nested function for psychoacoustic annoyance

end

