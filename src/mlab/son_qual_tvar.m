function [N, S, F, R] = son_qual_tvar(filepath, tnc, cal_val, cal_tnc,...
                                      cal_type, pctl)
% Returns time-varying sound quality metric values for an audio input file
% Calculates the time-varying loudness N according to ISO 532-1:2017 
% (Zwicker).
% Fluctuation strength and roughness are calculated based on Zwicker &
% Fastl, 1999. Psychoacoustics: Facts and Models (2nd edition).
% 
% Inputs
% ------
% filepath : string
%            the path to the input audio file
% tnc : vector or 2D array of floats, length 2
%       the truncation values (in seconds) to be applied to the start or
%       end of the input file - the first value in tnc applies to the start
%       and the second value to the end
% cal_val : double
%           the value of the calibration to be applied to the file
% cal_tnc : keyword string (default: 'pre')
%           specifies whether calibration of the input file should take
%           place pre- or post-truncation operation
%           'pre':apply calibration before file truncation
%           'post': apply calibration to the truncated file
% cal_type : keyword string (default: 'LAEQ')
%            specifies the type of calibration for the input file (ie, what
%            property cal_val represents)
% N_type : keyword string (default: '532-1')
%          '532-1':ISO 532-1:2017 / Zwicker loudness adapted for
%          time-varying signals
%           '532-2':ISO 532-2:2017 / Moore-Glasberg time-varying loudness
%           for binaural signals
% pctl : double or vector of integers 0-100 (default: 5)
%        the percentile to use in calculating an overall metric value from
%        the time-varying distribution
% 
% Returns
% -------
% N : structure containing percentile loudness, percentile value,
%       time-varying loudness, and specific time-varying loudness
%
% S : structure containing percentile sharpness, percentile value, and
%       time-varying sharpness
%
% R : structure containing percentile roughness, percentile value,
%       time-varying roughness, specific time-varying roughness, and roughness
%       modulation frequencies
%
% F : structure containing percentile fluctuation strength, percentile
%       value, time-varying fluctuation strength, specific time-varying
%       fluctuation strength, and fluctuation strength modulation frequencies
% Ownership and Quality Assurance
%
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%  
% Date created: 18/07/2023
% Date last modified: 19/07/2023
% MATLAB version: 2022b
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
    arguments
        filepath string
        tnc (1, 2) double {mustBePositive}
        cal_val (1, 1) double
        cal_tnc string {mustBeMember(cal_tnc, {'pre', 'post'})} = 'pre'
        cal_type string {mustBeMember(cal_type, {'LAEQ'})} = 'LAEQ'
        pctl (1, :) {mustBeNonnegative, mustBeInteger,...
                       mustBeLessThanOrEqual(pctl, 100)} = 5
    end

%% Input, truncation and calibration

[xt, fs] = audioread(filepath);

% signal calibration block (if pre-truncation)
if strcmpi(cal_tnc, 'pre')
    switch upper(cal_type)
        case 'LAEQ'
            weightFilt = weightingFilter('A-weighting', fs);
            xA_temp = weightFilt(xt);
            xA_rms = rms(xA_temp);
            cal_valPa = 2e-5*10^(cal_val/20);
            xn = xt.*(cal_valPa./xA_rms);
            clear xA_temp
    end
else
    xn = xt;
end

% signal truncation block
% start truncation
if tnc(1) ~= 0
    if tnc(1)*fs < size(xn, 1)
        tnc_s = floor(tnc(1)*fs);
        xn = xn(tnc_s+1:end, :);
    else
        error("Truncation value %i is larger than signal length %i", tnc(1),...
            size(xn, 1))
    end
end

% end truncation
if tnc(2) ~= 0
    if tnc(2)*fs < size(xn, 1)
        tnc_e = floor(tnc(2)*fs);
        xn = xn(1:end-tnc_e, :);
    else
        error("Truncation value %i is larger than signal length %i", tnc(2),...
            size(xn, 1))
    end
end

% signal calibration block (if post-truncation)
if strcmpi(cal_tnc, 'post')
    switch upper(cal_type) % calibrate to type
        case 'LAEQ'
            weightFilt = weightingFilter('A-weighting', fs);
            xA_temp = weightFilt(xn);
            xA_rms = rms(xA_temp);
            cal_valPa = 2e-5*10^(cal_val/20);
            x_Pa = xn.*(cal_valPa./xA_rms);
    end
else
    x_Pa = xn;
end

%% Signal processing

% loudness
[Nt, specNt, percN] = acousticLoudness(x_Pa, fs, TimeVarying=true,...
                                       Percentiles=pctl);

% sharpness
St = acousticSharpness(specNt, TimeVarying=true);
percS = prctile(St, pctl, 1);

% roughness
[Rt, specRt, fModR] = acousticRoughness(x_Pa, fs);
percR = prctile(Rt, pctl, 1);

% fluctuation strength
[Ft, specFt, fModF] = acousticFluctuation(specNt);
percF = prctile(Ft, pctl, 1);

% compile results into output data structures
% loudness
N = struct;
N.percN = percN;
N.pctl = pctl;
N.Nt = Nt;
N.specNt = specNt;

% sharpness
S = struct;
S.percS = percS;
S.pctl = pctl;
S.St = St;

% roughness
R = struct;
R.percR = percR;
R.pctl = pctl;
R.Rt = Rt;
R.specRt = specRt;
R.fModR = fModR;

% fluctation strength
F = struct;
F.percF = percF;
F.pctl = pctl;
F.Nt = Ft;
F.specFt = specFt;
F.fModF = fModF;

