function [f, X, Ben] = fft_aspec(xn, fs, Nfft, axisn, over, win, X_scale, twosided)
% [f, X, Ben] = fft_aspec(xn, fs, Nfft, axisn, over, win, X_side, X_scale)
% Returns FFT frequencies and magnitude auto spectra for input signal data 
% (NB: assumed to be stationary or quasi-stationary - transient signals may 
% require alternative processing).
% Uses Welch windowed averages to estimate the auto spectrum, which can be 
% used as the basis for other spectral scalings of the output:
% linear RMS spectrum = sqrt(autospectrum)
% linear peak spectrum = sqrt(2)*sqrt(autospectrum)
% power (auto) spectral density = autospectrum*(Ben/deltaFreq) where Ben is
% the normalised equivalent noise bandwidth for the window (see outputs
% below), and deltaFreq is the output spectral line spacing, f(2) - f(1).
% If one-sided, the returned spectral array X has a size like
% (Nfft/2 + odd, <xn.#signals>), where odd = 0 if Nfft is even, with the
% other dimension corresponding with the number of signals in xn (NB: the
% output X orientation matches the input orientation). If two-sided, the
% length of the spectrum/spectra will be Nfft.
%
% Inputs
% ------
% xn : vector or 2D array of floats
%      the input signals to be processed
% fs : integer
%      the sample frequency for the signal(s) - all signals must have same
%      fs
% axisn : integer 0-1 (default: 0, ie along rows)
%      the input axis over which to apply the FFT (time axis)
% Nfft : integer
%      the FFT block size - if this is smaller than the signal length,
%      ensemble averaging will be carried out using 'over' overlap
% over : decimal, >=0, < 1
%      the overlap proportion to use for windowed ensemble averaging
% win : keyword string (default: 'hann')
%      the window type to apply (supported: 'hann', 'hamming',
%      'flattopwin', 'rectwin'; aliases available for
%      'flattopwin':'flattop' and 'rectwin':'rect'|'boxcar')
% X_scale : keyword string (default: 'psd')
%      the type (scaling) of output (auto)spectrum. 
%      'psd':(auto) power spectral density (power spectrum/Hz);
%      'aspec':(auto) power spectrum;
%      'rms': linear root-mean-square amplitude spectrum;
%      'peak': linear peak amplitude spectrum
% twosided : Boolean
%            indicates the type of spectrum to output (one- or two-sided)
%
% Outputs
% -------
% f : vector
%      the spectral frequencies
% X : vector or 2D array
%      the real-valued autospectrum (or autospectra)
% Ben : double
%      the 'normalised equivalent noise bandwidth' for the window applied
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
% Date created: 03/01/2023
% Date last modified: 14/07/2024
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

%%  argument validation
    arguments (Input)
        xn (:, :) double
        fs (1, 1) {mustBePositive, mustBeInteger}
        Nfft (1, 1) {mustBePositive, mustBeInteger}
        axisn (1, 1) {mustBePositive, mustBeInteger,...
                       mustBeLessThanOrEqual(axisn, 2)} = 1
        over (1, 1) {mustBeNonnegative, mustBeLessThan(over, 1)} = 0.0
        win string {mustBeMember(win, {'hann', 'hamming', 'flattop',...
                                       'flattopwin', 'rectwin'...
                                       'rect', 'boxcar'})} = 'hann'
        X_scale string {mustBeMember(X_scale, {'psd', 'aspec', 'rms',...
                                               'peak'})} = 'psd'
        twosided {mustBeNumericOrLogical} = false
    end

    xnsize = size(xn); % dimensions of input signal array
    xnlength = xnsize(axisn); % signal length
    axes_i = [1, 2]; % axis indices
    xns = xnsize(axes_i(axes_i ~= axisn)); % number of signals in xn

    % check signal fits into block window
    if Nfft > xnlength
        error("FFT block length must not be longer than signal length.")
    end

%%  assign variable parameters
    sflag = 'periodic'; % asymmetry flag for windowing
    switch win
        case 'hann'
            fft_win = hann(Nfft, sflag);
        case 'hamming'
            fft_win = hamming(Nfft, sflag);
        case {'flattopwin', 'flattop'}
            fft_win = flattopwin(Nfft, sflag);
        case {'rectwin', 'rect', 'boxcar'}
            fft_win = rectwin(Nfft);
    end

    df = fs/Nfft;  % frequency interval for FFT spectrum
    Be = enbw(fft_win, fs);  % equivalent noise bandwidth for window
    Ben = Be/df;  % normalised equivalent noise bandwidth for window
    S_amp = Nfft/sum(fft_win);  % amplitude correction for window

    N = floor(Nfft/2);
    if mod(Nfft, 2) == 0
        odd = 0;
        f2 = df*[(0:N - 1), (-N:-1)];  % two-sided frequency vector for even Nfft
    else
        odd = 1;
        f2 = df*[(0:N), (-N:-1)];  % two-sided frequency vector for odd Nfft
    end
    
    if twosided
        f = f2;
    else
        f = f2(1:floor(Nfft/2) + odd);  % one-sided positive frequency vector
    end

    iiEnd = size(f);
    
%%  FFT processing
    % convert array to required form
    if axisn == 2
        xn2 = transpose(xn);
    else
        xn2 = xn;
    end
    axisn2 = 1;  % dummy axis variable

    % initialise unscaled autospectrum accumulation array
    XXn = zeros(Nfft, xns);
    block = 1;  % FFT processing block counter initialisation
    blstart = 1;  % block start index
    blend = Nfft;  % block end index

    while blend <= size(xn2, axisn2)
        % windowed block of xn
        xnw = repmat(fft_win, 1, xns).*xn2(blstart:blend, :);
        % unscaled autospectrum of windowed block
        XXn = XXn + conj(fft(xnw, Nfft, 1)).*fft(xnw, Nfft, 1);

        block = block + 1;  % advance block counter
        blstart = (block - 1)*Nfft*(1 - over) + 1;  % increment block start
        blend = blstart + Nfft - 1;  % increment block end
    end

    % total number of blocks
    nBlocks = block - 1;

    % average unscaled autospectrum accumulation array over blocks
    XX = XXn(1:length(f), :)./nBlocks;

    if ~twosided
        % discard negative frequencies (scaling to one-sided spectrum)
        XX(2:iiEnd, :) = 2*XX(2:iiEnd, :);
    else
        % shift negative frequencies
        f = fftshift(f);
        XX = fftshift(XX);

    % scale to autospectrum
    Axx = (S_amp/Nfft)^2*XX;

    switch X_scale
        case 'psd'
            X0 = Axx./Be;  % note: /df is incorporated into use of Be
                           % instead of Ben (= Be/df)
        case 'aspec'
            X0 = Axx;
        case 'rms'
            X0 = sqrt(Axx);
        case 'peak'
            X0 = sqrt(Axx)*sqrt(2);
    end
    
    % revert output array to input form
    if axisn == 2
        X = transpose(X0);  % non-conjungate tranpose
    else
        X = X0;
    end

end
