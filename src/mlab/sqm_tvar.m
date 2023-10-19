function [N, S, F, R] = sqm_tvar(p, fs, pcex)
% [N, S, F, R] = son_qual_tvar(p, fs, pcex)
% Returns time-varying sound quality metric values for an audio signal.
% Calculates the time-varying loudness N according to ISO 532-1:2017 
% (Zwicker).
% Fluctuation strength and roughness are calculated based on Zwicker &
% Fastl, 1999. Psychoacoustics: Facts and Models (2nd edition).
% 
% Inputs
% ------
% p : vector or 2D matrix of single mono or single stereo audio signals 
%     calibrated to sound pressure in Pascals (Pa)
%
% fs : integer
%       the sample frequency of the input signal(s)
%
% pcex : double or vector of integers 0-100 (default: 5)
%        the 'percent exceeded' value to use in calculating an overall metric
%        value from the time-varying distribution (note this is the percent
%        exceeded, and not the statistical percentile value of the
%        cumulative distribution function, which is 100 - pcex)
% 
% Returns
% -------
% N : structure containing percentile loudness, percentile value,
%       time-varying loudness, and specific time-varying loudness
%
% S : structure containing percentile sharpness, percentile value, and
%       time-varying sharpness
%
% F : structure containing percentile fluctuation strength, percentile
%       value, time-varying fluctuation strength, specific time-varying
%       fluctuation strength, and fluctuation strength modulation frequencies
%
% R : structure containing percentile roughness, percentile value,
%       time-varying roughness, specific time-varying roughness, and roughness
%       modulation frequencies
%
% Requirements
% ------------
% Audio Toolbox
%
% Ownership and Quality Assurance
%
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%  
% Date created: 18/07/2023
% Date last modified: 19/10/2023
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
        p (:, :) double {mustBeReal}
        fs (1, 1) double {mustBePositive, mustBeInteger}
        pcex (1, :) double {mustBeNonnegative, mustBeInteger,...
                            mustBeLessThanOrEqual(pcex, 100)} = 5
    end


%% Signal processing

% loudness
[Nt, specNt, percN] = acousticLoudness(p, fs, 1, TimeVarying=true,...
                                       Percentiles=pcex);

% sharpness
St = acousticSharpness(specNt, TimeVarying=true);
percS = prctile(St, 100 - pcex, 1);

% roughness
[Rt, specRt, fModR] = acousticRoughness(p, fs, 1);
percR = prctile(Rt, 100 - pcex, 1);

% fluctuation strength
[Ft, specFt, fModF] = acousticFluctuation(specNt);
percF = prctile(Ft, 100 - pcex, 1);

% compile results into output data structures
% loudness
N = struct;
N.percN = percN;
N.pctl = pcex;
N.Nt = Nt;
N.specNt = specNt;

% sharpness
S = struct;
S.percS = percS;
S.pctl = pcex;
S.St = St;

% roughness
R = struct;
R.percR = percR;
R.pctl = pcex;
R.Rt = Rt;
R.specRt = specRt;
R.fModR = fModR;

% fluctuation strength
F = struct;
F.percF = percF;
F.pctl = pcex;
F.Ft = Ft;
F.specFt = specFt;
F.fModF = fModF;

