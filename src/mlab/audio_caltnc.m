function [x_Pa, fs] = audio_caltnc(filepath, tnc, cal_val, cal_tnc, cal_type)
% [x_Pa, fs] = audio_caltnc(filepath, tnc, cal_val, cal_tnc, cal_type)
% Returns calibrated, truncated (optional) signal from audio input file.
% 
% Inputs
% ------
% filepath : string
%            the path to the input audio file
% tnc : vector or 2D matrix of floats, length 2
%       the truncation values (in seconds) to be applied to the start or
%       end of the input file - the first value in tnc applies to the start
%       and the second value to the end
% cal_val : double
%           the value of the calibration to be applied to the file
% cal_tnc : keyword string (default: 'pre')
%           specifies whether calibration of the input file should take
%           place pre- or post-truncation operation
%           'pre': apply calibration before file truncation
%           'post': apply calibration to the truncated file
% cal_type : keyword string (default: 'LAEQ')
%            specifies the type of calibration for the input file (ie, what
%            property cal_val represents)
%           'LAEQ': A-weighted energy time-averaged sound level
%           'LAE': A-weighted sound exposure level
%           'LPEAK': unweighted peak sound level
% 
% Returns
% -------
% x_Pa : vector or 2D matrix
%         the audio signal sample data calibrated to units of Pascals
%
% fs : positive integer
%      the sample frequency of the audio signal
%
% Requirements
% ------------
% Signal Processing Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%  
% Date created: 20/07/2023
% Date last modified: 20/07/2023
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
    arguments (Input)
        filepath string {mustBeFile}
        tnc (1, 2) {mustBeNonnegative}
        cal_val (1, 1) double
        cal_tnc string {mustBeMember(cal_tnc, {'pre', 'post'})} = 'pre'
        cal_type string {mustBeMember(cal_type, {'LAEQ', 'LAE',...
                                                  'LPEAK'})} = 'LAEQ'
    end

%% Input, truncation and calibration

[xt, fs] = audioread(filepath);

% signal calibration block (if pre-truncation)
if strcmpi(cal_tnc, 'pre')
    switch upper(cal_type)
        case 'LAEQ'
            % check fs is suitable for A filter
            if fs < 36000
                fs_re = 36000;
                up = fs_re/gcd(fs_re, fs);  % upsampling factor
                down = fs/gcd(fs_re, fs);  % downsampling factor
                xt_re = resample(xt, up, down);  % apply resampling
            else
                xt_re = xt;
                fs_re = fs;
            end
            weightFilt = weightingFilter('A-weighting', fs_re);
            xA_temp = weightFilt(xt_re);
            xA_rms = rms(xA_temp, 1);
            cal_valPa = 2e-5*10^(cal_val/20);
            xn = xt.*(cal_valPa./xA_rms);
            clear xA_temp
        case 'LAE'
            % check fs is suitable for A filter
            if fs < 36000
                fs_re = 36000;
                up = fs_re/gcd(fs_re, fs);  % upsampling factor
                down = fs/gcd(fs_re, fs);  % downsampling factor
                xt_re = resample(xt, up, down);  % apply resampling
            else
                xt_re = xt;
                fs_re = fs;
            end
            weightFilt = weightingFilter('A-weighting', fs_re);
            xA_temp = weightFilt(xt_re);
            xA_rss = sqrt(size(xA_temp, 1)/fs_re)*rms(xA_temp, 1);
            cal_valPa = 2e-5*10^(cal_val/20);
            xn = xt.*(cal_valPa./xA_rss);
            clear xA_temp
        case 'LPEAK'
            xpeak = max(abs(xt), [], 1);
            cal_valPa = 2e-5*10^(cal_val/20);
            xn = xt.*(cal_valPa./xpeak);
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
            size(xn, 1)/fs)
    end
end

% end truncation
if tnc(2) ~= 0
    if tnc(2)*fs < size(xn, 1)
        tnc_e = floor(tnc(2)*fs);
        xn = xn(1:end-tnc_e, :);
    else
        error("Truncation value %i is larger than signal length %i", tnc(2),...
            size(xn, 1)/fs)
    end
end

% signal calibration block (if post-truncation)
if strcmpi(cal_tnc, 'post')
    switch upper(cal_type) % calibrate to type
        case 'LAEQ'
            % check fs is suitable for A filter
            if fs < 36000
                fs_re = 36000;
                up = fs_re/gcd(fs_re, fs);  % upsampling factor
                down = fs/gcd(fs_re, fs);  % downsampling factor
                xn_re = resample(xn, up, down);  % apply resampling
            else
                xn_re = xn;
                fs_re = fs;
            end
            weightFilt = weightingFilter('A-weighting', fs_re);
            xA_temp = weightFilt(xn_re);
            xA_rms = rms(xA_temp, 1);
            cal_valPa = 2e-5*10^(cal_val/20);
            x_Pa = xn.*(cal_valPa./xA_rms);
            clear xA_temp
        case 'LAE'
            % check fs is suitable for A filter
            if fs < 36000
                fs_re = 36000;
                up = fs_re/gcd(fs_re, fs);  % upsampling factor
                down = fs/gcd(fs_re, fs);  % downsampling factor
                xn_re = resample(xn, up, down);  % apply resampling
            else
                xn_re = xn;
                fs_re = fs;
            end
            weightFilt = weightingFilter('A-weighting', fs_re);
            xA_temp = weightFilt(xn_re);
            xA_rss = sqrt(size(xA_temp, 1)/fs_re)*rms(xA_temp, 1);
            cal_valPa = 2e-5*10^(cal_val/20);
            x_Pa = xn.*(cal_valPa./xA_rss);
            clear xA_temp
        case 'LPEAK'
            xpeak = max(abs(xn), [], 1);
            cal_valPa = 2e-5*10^(cal_val/20);
            x_Pa = xn.*(cal_valPa./xpeak);
    end
else
    x_Pa = xn;
end
