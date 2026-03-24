function [signalFiltBank, fm] = iso5321_third_oct_filterbank(signal, sampleRate, axisN, fLim)
% ISO 532-1:2017 compliant 1/3-octave filter bank
%
% Return filtered signals.
%
% Inputs
% ------
% signal : vector or 2D matrix
%   Input signal as single or multichannel signals.
%
% sampleRate : integer
%   Sample rate (frequency) of the input signal(s).
%
% axisN : integer (1 or 2, default: 1)
%   Time axis along which to perform filtering.
%
% fLim : vector (default: [25, 12600])
%   Frequency limits for filter bank.
%
% Returns
% -------
% signalFiltBank : 2D or 3D matrix
%   Filtered signals arranged as [time, bands(, channels)].
% 
% fm : vector
%   Third-octave centre-frequencies for the output bands.
% 
% Assumptions
% -----------
%
% References
% ----------
% ISO 532-1:2017.
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
% Date created: 20/03/2026
% Date last modified: 23/03/2026
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
% Checked by:
% Date last checked:
%
    arguments (Input)
        signal (:, :) double {mustBeReal}
        sampleRate (1, 1) {mustBeInteger, mustBePositive} = 48e3
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        fLim (1, 2) double {mustBeInRange(fLim, 25, 12600)} = [25, 12600]
    end

% Ensure signal is orientated with column-wise time
if axisN == 2
    signal = signal.';
    % axisN = 1;
end

% Ensure sampling rate is full range
resampledRate = 48e3;
if sampleRate ~= resampledRate  % Resample signal
    up = resampledRate/gcd(resampledRate, sampleRate);  % upsampling factor
    down = sampleRate/gcd(resampledRate, sampleRate);  % downsampling factor
    signal = resample(signal, up, down);  % apply resampling
    % sampleRate = resampledRate;
end

[signalLen, signalChans] = size(signal);

% --- constants ---
nBandsTotal = 28;
nStages = 3;

idx = (0:nBandsTotal - 1);
fmAll = 1000*10.^((idx - 16)/10);

% --- select bands based on fLim ---
bandMask = (fmAll >= fLim(1)) & (fmAll <= fLim(2));
bandIdx = find(bandMask);

fm = fmAll(bandIdx);
nBands = numel(bandIdx);

% Reference coefficients
ref = [1, 2, 1,  1, -2, 1;
       1, 0, -1,  1, -2, 1;
       1, -2, 1,  1, -2, 1];

% Load the difference table
diffTab = getThirdOctaveDiffTable();
gainTab = getThirdOctaveGainTable();

% --- precompute coefficients ---
coeffAll = zeros(nBandsTotal, nStages, 6);
for b = 1:nBandsTotal
    for s = 1:nStages
        coeffAll(b, s, :) = ref(s, :) - squeeze(diffTab(b, s, :))';
    end
end

% Initialise output
signalFiltBank = zeros(signalLen, nBands, signalChans);

% --- reusable buffers ---
stageInput  = zeros(signalLen, signalChans);
stageOutput = zeros(signalLen, signalChans);

% Loop over bands to apply filter
for bb = 1:nBands

    band = bandIdx(bb);

    stageInput(:, :) = signal;

    for stage = 1:nStages

        % Get coefficient values
        c = squeeze(coeffAll(band, stage, :)).';
        G = gainTab(band, stage);

        b0 = c(1); b1 = c(2); b2 = c(3);
        a1 = c(5); a2 = c(6);

        % Reset states
        w1 = zeros(1, signalChans);
        w2 = zeros(1, signalChans);

        % Filter
        for n = 1:signalLen
            xn = stageInput(n, :);

            w0 = xn*G - a1*w1 - a2*w2;
            stageOutput(n, :) = b0*w0 + b1*w1 + b2*w2;

            w2 = w1;
            w1 = w0;
        end

        % Buffer swap
        tmp = stageInput;
        stageInput = stageOutput;
        stageOutput = tmp;
    end

    signalFiltBank(:, bb, :) = stageInput;
end


function gainTab = getThirdOctaveGainTable()
    
    gainTab = [ ...
        4.30764e-011 1 1;
        8.59340e-011 1 1;
        1.71424e-010 1 1;
        3.41944e-010 1 1;
        6.82035e-010 1 1;
        1.36026e-009 1 1;
        2.71261e-009 1 1;
        5.40870e-009 1 1;
        1.07826e-008 1 1;
        2.14910e-008 1 1;
        4.28228e-008 1 1;
        8.54316e-008 1 1;
        1.70009e-007 1 1;
        3.38215e-007 1 1;
        6.71990e-007 1 1;
        1.33531e-006 1 1;
        2.65172e-006 1 1;
        5.25477e-006 1 1;
        1.03780e-005 1 1;
        2.04870e-005 1 1;
        4.05198e-005 1 1;
        7.97914e-005 1 1;
        1.56511e-004 1 1;
        3.04954e-004 1 1;
        5.99157e-004 1 1;
        1.16544e-003 1 1;
        2.27488e-003 1 1;
        3.91006e-003 1 1];

end

function diffTab = getThirdOctaveDiffTable()

    diffTab = zeros(28,3,6);
    
    diffTab(1,:,:) = [ ...
        0 0 0 0 -6.70260e-004  6.59453e-004;
        0 0 0 0 -3.75071e-004  3.61926e-004;
        0 0 0 0 -3.06523e-004  2.97634e-004];
    
    diffTab(2,:,:) = [ ...
        0 0 0 0 -8.47258e-004  8.30131e-004;
        0 0 0 0 -4.76448e-004  4.55616e-004;
        0 0 0 0 -3.88773e-004  3.74685e-004];
    
    diffTab(3,:,:) = [ ...
        0 0 0 0 -1.07210e-003  1.04496e-003;
        0 0 0 0 -6.06567e-004  5.73553e-004;
        0 0 0 0 -4.94004e-004  4.71677e-004];
    
    diffTab(4,:,:) = [ ...
        0 0 0 0 -1.35836e-003  1.31535e-003;
        0 0 0 0 -7.74327e-004  7.22007e-004;
        0 0 0 0 -6.29154e-004  5.93771e-004];
    
    diffTab(5,:,:) = [ ...
        0 0 0 0 -1.72380e-003  1.65564e-003;
        0 0 0 0 -9.91780e-004  9.08866e-004;
        0 0 0 0 -8.03529e-004  7.47455e-004];
    
    diffTab(6,:,:) = [ ...
        0 0 0 0 -2.19188e-003  2.08388e-003;
        0 0 0 0 -1.27545e-003  1.14406e-003;
        0 0 0 0 -1.02976e-003  9.40900e-004];
    
    diffTab(7,:,:) = [ ...
        0 0 0 0 -2.79386e-003  2.62274e-003;
        0 0 0 0 -1.64828e-003  1.44006e-003;
        0 0 0 0 -1.32520e-003  1.18438e-003];
    
    diffTab(8,:,:) = [ ...
        0 0 0 0 -3.57182e-003  3.30071e-003;
        0 0 0 0 -2.14252e-003  1.81258e-003;
        0 0 0 0 -1.71397e-003  1.49082e-003];
    
    diffTab(9,:,:) = [ ...
        0 0 0 0 -4.58305e-003  4.15355e-003;
        0 0 0 0 -2.80413e-003  2.28135e-003;
        0 0 0 0 -2.23006e-003  1.87646e-003];
    
    diffTab(10,:,:) = [ ...
        0 0 0 0 -5.90655e-003  5.22622e-003;
        0 0 0 0 -3.69947e-003  2.87118e-003;
        0 0 0 0 -2.92205e-003  2.36178e-003];
    
    diffTab(11,:,:) = [ ...
        0 0 0 0 -7.65243e-003  6.57493e-003;
        0 0 0 0 -4.92540e-003  3.61318e-003;
        0 0 0 0 -3.86007e-003  2.97240e-003];
    
    diffTab(12,:,:) = [ ...
        0 0 0 0 -1.00023e-002  8.29610e-003;
        0 0 0 0 -6.63788e-003  4.55999e-003;
        0 0 0 0 -5.15982e-003  3.75306e-003];
    
    diffTab(13,:,:) = [ ...
        0 0 0 0 -1.31230e-002  1.04220e-002;
        0 0 0 0 -9.02274e-003  5.73132e-003;
        0 0 0 0 -6.94543e-003  4.71734e-003];
    
    diffTab(14,:,:) = [ ...
        0 0 0 0 -1.73693e-002  1.30947e-002;
        0 0 0 0 -1.24176e-002  7.20526e-003;
        0 0 0 0 -9.46002e-003  5.93145e-003];
    
    diffTab(15,:,:) = [ ...
        0 0 0 0 -2.31934e-002  1.64308e-002;
        0 0 0 0 -1.73009e-002  9.04761e-003;
        0 0 0 0 -1.30358e-002  7.44926e-003];
    
    diffTab(16,:,:) = [ ...
        0 0 0 0 -3.13292e-002  2.06370e-002;
        0 0 0 0 -2.44342e-002  1.13731e-002;
        0 0 0 0 -1.82108e-002  9.36778e-003];
    
    diffTab(17,:,:) = [ ...
        0 0 0 0 -4.28261e-002  2.59325e-002;
        0 0 0 0 -3.49619e-002  1.43046e-002;
        0 0 0 0 -2.57855e-002  1.17912e-002];
    
    diffTab(18,:,:) = [ ...
        0 0 0 0 -5.91733e-002  3.25054e-002;
        0 0 0 0 -5.06072e-002  1.79513e-002;
        0 0 0 0 -3.69401e-002  1.48094e-002];
    
    diffTab(19,:,:) = [ ...
        0 0 0 0 -8.26348e-002  4.05894e-002;
        0 0 0 0 -7.40348e-002  2.24476e-002;
        0 0 0 0 -5.34977e-002  1.85371e-002];
    
    diffTab(20,:,:) = [ ...
        0 0 0 0 -1.17018e-001  5.08116e-002;
        0 0 0 0 -1.09516e-001  2.81387e-002;
        0 0 0 0 -7.85097e-002  2.32872e-002];
    
    diffTab(21,:,:) = [ ...
        0 0 0 0 -1.67714e-001  6.37872e-002;
        0 0 0 0 -1.63378e-001  3.53729e-002;
        0 0 0 0 -1.16419e-001  2.93723e-002];
    
    diffTab(22,:,:) = [ ...
        0 0 0 0 -2.42528e-001  7.98576e-002;
        0 0 0 0 -2.45161e-001  4.43370e-002;
        0 0 0 0 -1.73972e-001  3.70015e-002];
    
    diffTab(23,:,:) = [ ...
        0 0 0 0 -3.53142e-001  9.96330e-002;
        0 0 0 0 -3.69163e-001  5.53535e-002;
        0 0 0 0 -2.61399e-001  4.65428e-002];
    
    diffTab(24,:,:) = [ ...
        0 0 0 0 -5.16316e-001  1.24177e-001;
        0 0 0 0 -5.55473e-001  6.89403e-002;
        0 0 0 0 -3.93998e-001  5.86715e-002];
    
    diffTab(25,:,:) = [ ...
        0 0 0 0 -7.56635e-001  1.55023e-001;
        0 0 0 0 -8.34281e-001  8.58123e-002;
        0 0 0 0 -5.94547e-001  7.43960e-002];
    
    diffTab(26,:,:) = [ ...
        0 0 0 0 -1.10165e+000  1.91713e-001;
        0 0 0 0 -1.23939e+000  1.05243e-001;
        0 0 0 0 -8.91666e-001  9.40354e-002];
    
    diffTab(27,:,:) = [ ...
        0 0 0 0 -1.58477e+000  2.39049e-001;
        0 0 0 0 -1.80505e+000  1.28794e-001;
        0 0 0 0 -1.32500e+000  1.21333e-001];
    
    diffTab(28,:,:) = [ ...
        0 0 0 0 -2.50630e+000  1.42308e-001;
        0 0 0 0 -2.19464e+000  2.76470e-001;
        0 0 0 0 -1.90231e+000  1.47304e-001];

end

end