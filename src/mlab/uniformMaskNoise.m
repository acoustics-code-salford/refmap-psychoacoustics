function [x, h] = uniformMaskNoise(duration, fs, levelRMS)
% uniformMaskNoise  Uniform masking noise (Fastl & Zwicker) = white noise shaped by the
%      attenuation curve below (0 dB to ~500 Hz, ~15 dB at 20 kHz).
%   [x,h] = umn(duration_s, fs, levelRMS)   fs >= 40 kHz; levelRMS optional.

    if nargin < 2, fs = 48000; end

    % Attenuation curve: [Hz, dB]
    curve = [  20    0
              500    0
              700  1
             1000  2
             2000  4
             5000  8
            10000 12
            20000 15];

    numtaps = 4096;                                  % fir2 order (filter has order+1 taps)
    f = linspace(0, fs/2, 2048);                     % Hz
    fc = min(max(f, curve(1,1)), curve(end,1));      % hold end values outside the curve
    attenDB = interp1(log10(curve(:,1)), curve(:,2), log10(fc), 'pchip');
    gain = 10.^(-attenDB/20);                        % attenuation -> linear gain
    h = fir2(numtaps, f/(fs/2), gain);               % linear-phase FIR

    n = round(duration*fs);
    white = randn(n + numtaps, 1);
    y = filter(h, 1, white);
    x = y(numtaps+1 : end);                          % discard filter start-up

    paRMS = 2e-5*(levelRMS/20);

    if nargin >= 3 && ~isempty(paRMS)
        x = x * paRMS / rms(x);
    end
end