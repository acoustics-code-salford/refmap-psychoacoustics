function [f, X, Ben] = rfft_aspec(xr, fs, Nfft, axisn, over, win, X_scale)
% [f, X, Ben] = rfft_aspec(xr, fs, Nfft, axisn, over, win, X_scale)
% Returns one-sided FFT frequencies and magnitude auto spectra for real-valued
% input signal data (NB: assumed to be stationary or quasi-stationary
% - transient signals may require alternative processing).
% Uses Welch windowed averages to estimate the auto spectrum, which can be 
% used as the basis for other spectral scalings of the output:
% linear RMS spectrum = sqrt(autospectrum)
% linear peak spectrum = sqrt(2)*sqrt(autospectrum)
% power (auto) spectral density = autospectrum*(Ben/deltaFreq) where Ben is
% the normalised equivalent noise bandwidth for the window (see outputs
% below), and deltaFreq is the output spectral line spacing.
% The returned spectral array X has a size like (Nfft/2 + 1, <xr.#signals>),
% with the other dimension corresponding with the number of signals in xr
% (NB: the output X orientation matches the input orientation).
%
% Inputs
% ------
% xr : vector or 2D array of real-valued floats
%      the input signals to be processed
% fs : integer
%      the sample frequency for the signal(s) - all signals must have same
%      fs
% axisn : integer 0-1 (default: 0, ie along rows)
%         the input axis over which to apply the FFT (time axis)
% Nfft : integer
%        the FFT block size - if this is smaller than the signal length,
%        ensemble averaging will be carried out using 'over' overlap
% over : float, >=0, < 1
%       the overlap proportion to use for windowed ensemble averaging
% win : keyword string (default: 'hann')
%       the window type to apply (supported: 'hann', 'hamming',
%       'flattopwin', 'rectwin'; aliases available for
%       'flattopwin':'flattop' and 'rectwin':'rect'|'boxcar')
% X_scale : keyword string (default: 'psd')
%           the type (scaling) of output (auto)spectrum. 
%           'psd':(auto) power spectral density (power spectrum/Hz);
%           'aspec':(auto) power spectrum;
%           'rms': linear root-mean-square amplitude spectrum;
%           'peak': linear peak amplitude spectrum
%
% Outputs
% -------
% f : vector
%      the one-sided (non-negative) spectral frequencies
% X : vector or 2D array
%      the real-valued, one-sided autospectrum (or autospectra)
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
        xr (:, :) {mustBeReal}
        fs (1, 1) {mustBePositive, mustBeInteger}
        Nfft (1, 1) {mustBePositive, mustBeInteger}
        axisn (1, 1) {mustBePositive, mustBeInteger,...
                       mustBeLessThanOrEqual(axisn, 2)} = 1
        over (1, 1) {mustBeNonnegative, mustBeLessThan(over, 1)} = 0.0
        win string {mustBeMember(win, {'hann', 'hamming', 'flattop', ...
                                       'flattopwin', 'rectwin'...
                                       'rect', 'boxcar'})} = 'hann'
        X_scale string {mustBeMember(X_scale, {'psd', 'aspec', 'rms',...
                                               'peak'})} = 'psd'
    end

    xrsize = size(xr); % dimensions of input signal array
    xrlength = xrsize(axisn); % signal length
    axes_i = [1, 2]; % axis indices
    xrs = xrsize(axes_i(axes_i ~= axisn)); % number of signals in xr

    % check signal fits into block window
    if Nfft > xrlength
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

    f = df*(0:floor(Nfft/2));  % one-sided frequency vector
                                % (including k=N/2 for even inputs)

%%  FFT processing
    % convert array to required form
    if axisn == 2
        xr2 = transpose(xr);
    else
        xr2 = xr;
    end
    axisn2 = 1;  % dummy axis variable

    % initialise unscaled autospectrum accumulation array
    XXn = zeros(Nfft, xrs);
    block = 1;  % FFT processing block counter initialisation
    blstart = 1;  % block start index
    blend = Nfft;  % block end index

    while blend <= size(xr2, axisn2)
        % windowed block of xr
        xrw = repmat(fft_win, 1, xrs).*xr2(blstart:blend, :);
        % unscaled autospectrum of windowed block
        XXn = XXn + conj(fft(xrw, Nfft, 1)).*fft(xrw, Nfft, 1);

        block = block + 1;  % advance block counter
        blstart = (block - 1)*Nfft*(1 - over) + 1;  % increment block start
        blend = blstart + Nfft - 1;  % increment block end
    end

    % total number of blocks
    nBlocks = block - 1;

    % average unscaled autospectrum accumulation array over blocks and
    % discard negative frequencies (scaling to one-sided spectrum)
    % (NB: the k=N/2 line is retained, leading to a (Nfft/2)+1 length spectrum)
    XX = XXn(1:length(f), :)./nBlocks;
    XX(2:end, :) = 2*XX(2:end, :);
    
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
