function [Tt, ft_ton, Ttz, fz_ton, T, Tz] = acousticHMSTonality(p, r_s, axisn, outplot, ecma)
% [Tt, ft_ton, Ttz, fz_ton, T, Tz] = acousticHMTonality(p, r_s, axisn, outplot)
% Returns tonality values according to ECMA-418-2:2022 (using the Hearing
% Model according to Sottek) for an input calibrated single mono or single
% stereo audio (sound pressure) time-series signal, p.
%
% Inputs
% ------
% p : vector or 2D matrix
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% r_s : integer
%      the sample rate (frequency) of the input signal(s)
%
% axisn : integer (1 or 2, default: 1)
%         the time axis along which to calculate the tonality
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% ecma : Boolean true/false (default: true)
%        flag indicating whether to maintain strict standard adherence to
%        ECMA-418-2:2022 Equation 40, or otherwise to use an alternative
%        that provides closer time-alignment of the time-dependent tonality
%        with the original signal
% 
% Returns
% -------
% Tt : vector or 2D matrix
%      time-dependent total tonality values for each input signal channel
%
% ft_ton : vector or 2D matrix
%          time-dependent frequencies of the dominant tonal components
%          corresponding with the time-dependent total tonality values for
%          each input signal channel
%
% Ttz : vector or 2D matrix
%       time-dependent specific tonality values in each half-critical band
%       rate scale width for each input signal channel
%
% ftz_ton : vector or 2D matrix
%           time-dependent frequencies of the dominant tonal components
%           corresponding with each of the time-dependent specific tonality
%           values in each half-critical band rate scale width for each
%           input signal channel
%
% T : number or vector
%     total tonality value for each channel in the input signal
% 
% Tz : vector or 2D matrix
%      specific tonality values for each input signal channel
%
% Assumptions
% -----------
% None
%
% Requirements
% ------------
% Signal Processing Toolbox
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (matt@anv.co.uk)
% Institution: University of Salford / ANV Measurement Systems
%
% Date created: 07/08/2023
% Date last modified: 19/09/2023
% MATLAB version: 2022b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
%
% This code was developed from an original file 'SottekTonality.m' authored
% by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
% The original code has been reused and updated here with permission.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        p (:, :) double {mustBeReal}
        r_s (1, 1) double {mustBePositive, mustBeInteger}
        axisn (1, 1) {mustBeInteger, mustBeInRange(axisn, 1, 2)} = 1
        outplot {mustBeNumericOrLogical} = false
        ecma {mustBeNumericOrLogical} = true
    end

% Orient input matrix
if axisn == 2
    p = p.';
end

% Check the length of the input data
if size(p, 1) <  14336
    error('Error: Input signal is too short to calculate tonality')
    return
end

% Check the channel number of the input data
if size(p, 2) > 2
    error('Error: Input signal comprises more than two channels')
    return
end

%% Define constants

r_sre = 48e3;  % Signal sample rate prescribed to be 48kHz, Section 5.1.1 ECMA-418-2:2022
df0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2022
c = 0.1618;  % defined in Section 5.1.4.1 ECMA-418-2:2022

z = 0.5:0.5:26.5;  % half-critical band rate scale
Fz = (df0/c)*sinh(c*z);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2022
dfz = sqrt(df0^2 + (c*Fz).^2);  % Section 5.1.4.1 Equation 10 ECMA-418-2:2022

c_N = 0.0211668;  % Calibration factor from Section 5.1.8 Equation 23 ECMA-418-2:2022
c_x = 1.0005;  % Adjustment to calibration factors Footnotes 9, 22 & 33 ECMA-418-2:2022
a = 1.5;  % Constant from Section 5.1.8 Equation 23 ECMA-418-2:2022

% Values from Section 5.1.8 Table 2 ECMA-418-2:2022
p_t = 2e-5*10.^((15:10:85)/20).';
v = [1, 0.6602, 0.0864, 0.6384, 0.0328, 0.4068, 0.2082, 0.3994, 0.6434];

% Loudness threshold in quiet Section 5.1.9 Table 3 ECMA-418-2:2022
LTQz = [0.3310, 0.1625, 0.1051, 0.0757, 0.0576, 0.0453, 0.0365, 0.0298,...
        0.0247, 0.0207, 0.0176, 0.0151, 0.0131, 0.0115, 0.0103, 0.0093,...
        0.0086, 0.0081, 0.0077, 0.0074, 0.0073, 0.0072, 0.0071, 0.0072,...
        0.0073, 0.0074, 0.0076, 0.0079, 0.0082, 0.0086, 0.0092, 0.0100,...
        0.0109, 0.0122, 0.0138, 0.0157, 0.0172, 0.0180, 0.0180, 0.0177,...
        0.0176, 0.0177, 0.0182, 0.0190, 0.0202, 0.0217, 0.0237, 0.0263,...
        0.0296, 0.0339, 0.0398, 0.0485, 0.0622];

% Block and hop sizes Section 6.2.2 Table 4 ECMA-418-2:2022
overlap = 0.75;  % block overlap proportion
% block sizes
sz_b = [8192*ones(1, 3), 4096*ones(1, 13), 2048*ones(1, 9), 1024*ones(1, 28)];
sz_h = (1 - overlap)*sz_b;  % hop sizes (section 5.1.2 footnote 3 ECMA 418-2:2022)

r_sd = r_sre/min(sz_h);  % Output sample rate based on hop sizes - Resampling to common time basis Section 6.2.6 ECMA-418-2:2022

% Number of bands that need averaging. Section 6.2.3 Table 5 ECMA-418-2:2022
NB = [0, 1, 2*ones(1,14), ones(1,9), zeros(1,28);...
      1, 1, 2*ones(1,14), ones(1,9), zeros(1,28)];

% Critical band interpolation factors from Section 6.2.6 Table 6 ECMA-418-2:2022
i_interp = sz_b/min(sz_b);

% Noise reduction constants from Section 6.2.7 Table 7 ECMA-418-2:2022
alpha = 20;
beta = 0.07;

% Sigmoid function factor parameters Section 6.2.7 Table 8 ECMA-418-2:2022
csz_b = [18.21*ones(1, 3), 12.14*ones(1, 13), 417.54*ones(1, 9),...
         962.68*ones(1, 28)]; 
dsz_b = [0.36*ones(1, 3), 0.36*ones(1, 13), 0.71*ones(1, 9),...
         0.69*ones(1, 28)]; 

% Scaling factor constants from Section 6.2.8 Table 9 ECMA-418-2:2022
A = 35;
B = 0.003;
c_T = 2.8785151;  % calibration factor in Section 6.2.8 Equation 51 ECMA-418-2:2022

%% Signal processing

% Input pre-processing
% --------------------
if r_s ~= r_sre  % Resample signal
    up = r_sre/gcd(r_sre, r_s);  % upsampling factor
    down = r_s/gcd(r_sre, r_s);  % downsampling factor
    p_re = resample(p, up, down);  % apply resampling
else  % don't resample
    p_re = p;
end

% Fade in weighting function Section 5.1.2 ECMA-418-2:2022
w_fade = repmat(transpose(0.5 - 0.5*cos(pi*(0:239)/240)), 1, size(p_re, 2));
% Apply fade in
p_rew = [w_fade.*p_re(1:240, :);
         p_re(241:end, :)];

% Zero-padding Section 5.1.2 ECMA-418-2:2022
n_zeross = max(sz_b);  % start padding
n_samples = size(p_re, 1);
n_new = max(sz_h)*(ceil((n_samples + max(sz_h) + n_zeross)/max(sz_h)) - 1);
n_zerose = n_new - n_samples;  % end padding
% Apply padding
pn = [zeros(n_zeross, size(p_rew, 2));
      p_rew;
      zeros(n_zerose, size(p_rew, 2))];

% Apply outer & middle ear filter bank
% ------------------------------------
%
% Filter coefficients from Section 5.1.3.2 Table 1 ECMA-418-2:2022
% b_0k = [1.015896, 0.958943, 0.961372, 2.225804, 0.471735, 0.115267, 0.988029,...
%         1.952238];
% b_1k = [-1.925299, -1.806088, -1.763632, -1.43465, -0.366092, 0.0, -1.912434,...
%         0.16232]; 
% b_2k = [0.922118, 0.876439, 0.821788, -0.498204, 0.244145, -0.115267,...
%         0.926132, -0.667994];
% a_0k = ones(size(b_0k));
% a_1k = [-1.925299, -1.806088, -1.763632, -1.43465, -0.366092, -1.796003,...
%         -1.912434, 0.16232];
% a_2k = [0.938014, 0.835382, 0.78316, 0.727599, -0.28412, 0.805838, 0.914161,...
%         0.284244]; 
%
% Accurate coefficient values
b_0k = [1.01589602025559, 0.958943219304445, 0.961371976333197,...
        2.22580350360974, 0.471735128494163, 0.115267139824401,...
        0.988029297230954, 1.95223768730136];
b_1k = [-1.92529887777608, -1.80608801184949, -1.76363215433825,...
        -1.43465048479216, -0.366091796830044, 0.0, -1.91243380293387,...
        0.162319983017519];
b_2k = [0.922118060364679, 0.876438777856084, 0.821787991845146,...
        -0.498204282194628, 0.244144703885020, -0.115267139824401,...
        0.926131550180785, -0.667994113035186];
a_0k = ones(size(b_0k));
a_1k = [-1.92529887777608, -1.80608801184949, -1.76363215433825,...
        -1.43465048479216, -0.366091796830044, -1.79600256669201,...
        -1.91243380293387, 0.162319983017519];
a_2k = [0.938014080620272, 0.835381997160530, 0.783159968178343,...
        0.727599221415107, -0.284120167620817, 0.805837815618546,...
        0.914160847411739, 0.284243574266175]; 

sos = [b_0k.', b_1k.', b_2k.', a_0k.', a_1k.', a_2k.'];

% Section 5.1.3.2 ECMA-418-2:2022 Outer and middle/inner ear filtering
pn_om = sosfilt(sos, pn, 1);

% Check channels
nchans = size(pn_om, 2);
if nchans > 1
    chans = ["Stereo left";
             "Stereo right"];
else
    chans = "Mono";
end

n_steps = 115;  % approximate number of calculation steps
% Loop through channels in file
% -----------------------------
for chan = size(pn_om, 2):-1:1
    w = waitbar(0, "Initialising...");
    step_i = 1;
    % Apply auditory filter bank
    % --------------------------
    waitbar(step_i/n_steps, w, 'Applying auditory filters...');
    step_i = step_i + 1;
    % Filter equalised signal using 53 1/2bark ERB filters according to 
    % Section 5.1.4.2 ECMA-418-2:2022

    k = 5;  % filter order = 5, footnote 5 ECMA-418-2:2022
    e_i = [0, 1, 11, 11, 1];  % filter coefficients for Section 5.1.4.2 Equation 15 ECMA-418-2:2022
    
    for zz = 53:-1:1
        % Section 5.1.4.1 Equation 8 ECMA-418-2:2022
        tau = (1/(2^(2*k - 1))).*nchoosek(2*k - 2, k - 1).*(1./dfz(zz));
        
        d = exp(-1./(r_sre.*tau)); % Section 5.1.4.1 ECMA-418-2:2022
        
        % Band-pass modifier Section 5.1.4.2 Equation 16/17 ECMA-418-2:2022
        bp = exp((1i.*2.*pi.*Fz(zz).*(0:k+1))./r_sre);
        
        % Feed-backward coefficients, Section 5.1.4.2 Equation 14 ECMA-418-2:2022
        m = 1:k;
        a_m = ([1, ((-d).^m).*arrayfun(@(m_) nchoosek(k, m_), m)]).*bp(1:k+1);
    
        % Feed-forward coefficients, Section 5.1.4.2 Equation 15 ECMA-418-2:2022
        m = 0:k-1;
        i = 1:k-1;
        b_m = ((((1-d).^k)./sum(e_i(i+1).*(d.^i))).*(d.^m).*e_i).*bp(1:k);
    
        % Recursive filter Section 5.1.4.2 Equation 13 ECMA-418-2:2022
        % Note, the results are complex so 2x the real-valued band-pass signal
        % is required.
        pn_omz(:, zz) = 2*real(filter(b_m, a_m, pn_om(:, chan)));

    end
    
    % Half Wave Rectification
    % -----------------------
    % Section 5.1.6 Equation 16 ECMA-418-2:2020
    pn_romz = pn_omz;
    pn_romz(pn_romz <= 0) = 0;

    % Autocorrelation function analysis
    % ---------------------------------
    % Duplicate Banded Data for ACF
    % Averaging occurs over neighbouring bands, to do this the segmentation
    % needs to be duplicated for neigbouring bands.  '_temp' has been added
    % to variables to indicate that the vectors/matrices have been modified
    % for duplicated neigbouring bands.
    
    pn_romz_temp = [pn_romz(:, 1:5), pn_romz(:, 2:18), pn_romz(:, 16:26), pn_romz(:, 26:53)];
    sz_b_temp = [8192*ones(1, 5) 4096*ones(1, 17) 2048*ones(1, 11) 1024*ones(1, 28)];
    sz_h_temp = (1 - overlap)*sz_b_temp;
    LTQz_temp = [LTQz(1:5), LTQz(2:18), LTQz(16:26), LTQz(26:53)];
    NBi_temp = [1, 1, 1, 6:18, 23:31, 34:61;
                2, 3, 5, 10:22, 25:33, 34:61];  % these are the (duplicated) indices corresponding with the NB bands around each z band
    
    for zz = 61:-1:1
    
        waitbar(((62 - zz) + step_i)/n_steps, w, strcat("Applying ACF in 61 bands, ",...
            num2str(zz), ' to go...'));
        
        % Section 5.1.5 Equation 19 ECMA-418-2:2022
        i_start = sz_b(1) - sz_b_temp(zz) + 1;

        % Truncate the signal to start from i_start and to end at an index
        % corresponding with the truncated signal length that will fill an
        % integer number of overlapped blocks
        pn_romzt = pn_romz_temp(i_start:end, zz);
        n_blocks = floor((size(pn_romzt, 1) - overlap*sz_b_temp(zz))/(sz_b_temp(zz) - overlap*sz_b_temp(zz)));
        i_end = n_blocks*sz_b_temp(zz)*(1 - overlap) + overlap*sz_b_temp(zz);
        pn_romzt = pn_romzt(1:i_end);

        % Arrange the signal into overlapped blocks - each block reads
        % along first axis, and each column is the succeeding overlapped
        % block. 3 columns of zeros are appended to the left side of the
        % matrix and the column shifted copies of this matrix are
        % concatenated. The first 6 columns are then discarded as these all
        % contain zeros from the appended zero columns.
        pn_lz = [zeros(sz_h_temp(zz), 3), reshape(pn_romzt, sz_h_temp(zz), [])];

        pn_lz = cat(1, circshift(pn_lz, 3, 2), circshift(pn_lz, 2, 2),...
                  circshift(pn_lz, 1, 2), circshift(pn_lz, 0, 2));

        pn_lz = pn_lz(:, 7:end);
  
        % Apply ACF
        % ACF implementation using DFT
        % Section 6.2.2 Equations 27 & 28 ECMA-418-2:2022
        phim_ulz = ifft(abs(fft(pn_lz, 2*sz_b_temp(zz), 1)).^2, 2*sz_b_temp(zz), 1);
        % Section 6.2.2 Equation 29 ECMA-418-2:2022
        denom = sqrt(cumsum(pn_lz.^2, 1, 'reverse').*flipud(cumsum(pn_lz.^2))) + 1e-12; 
        phim_lz = phim_ulz(1:sz_b_temp(zz), :)./denom;  % note that the block length is used here, rather than the 2*s_b, for compatability with the remaining code - beyond 0.75*s_b is assigned (unused) zeros in the next line
        phim_lz((0.75*sz_b_temp(zz) + 1):sz_b_temp(zz), :) = 0;

        % Transformation into Loudness
        % ----------------------------
        % Section 5.1.7 Equation 22 ECMA-418-2:2022
        plz{zz} = sqrt((2/sz_b_temp(zz))*sum(pn_lz.^2, 1));
        % Section 5.1.8 Equations 23 & 24 ECMA-418-2:2022
        Nz = c_N*c_x*(plz{zz}/20e-6).*prod((1 + (plz{zz}./p_t).^a).^((diff(v)/a)'));
        % Section 5.1.9 Equation 25 ECMA-418-2:2022
        Nz_basis_z = Nz - LTQz_temp(zz);
        Nz_basis_z((Nz-LTQz_temp(zz)) < 0) = 0;
        Nz_basis{zz} = Nz_basis_z;

        % Section 6.2.2 Equation 30 ECMA-418-2:202
        phim_lz_temp{zz} = Nz_basis_z.*phim_lz;

    end
    
    step_i = step_i + 62;  % increment calculation step for waitbar

    % Average the ACF over nB bands - Section 6.2.3 ECMA-418-2:2022        
    for zz = 53:-1:1  % Loop through 53 critical band filtered signals
        waitbar(((54 - zz) + step_i)/n_steps, w, strcat("Calculating tonality in 53 bands, ",...
            num2str(zz), " to go..."));
        
        NBZ = NB(1, zz) + NB(2, zz) + 1; % Total number of bands to average over
        
        % Averaging of frequency bands
        phim_lz_avg = mean(reshape(cell2mat(phim_lz_temp(NBi_temp(1, zz):NBi_temp(2, zz))), sz_b(zz), [], NBZ), 3);

        % Average the ACF over adjacent time blocks
        if zz <= 16 
            phim_lz_avg = movmean(phim_lz_avg, 3, 2, 'omitnan', 'EndPoints', 'fill');
        end
        
        % Application of ACF lag window Section 6.2.4 ECMA-418-2:2022
        tauz_start = max(0.5/dfz(zz), 2e-3);  % Equation 31 ECMA-418-2:2022
        tauz_end = max(4/dfz(zz), tauz_start + 1e-3);  % Equation 32 ECMA-418-2:2022
        % Equations 33 & 34 ECMA-418-2:2022
        mz_start = ceil(tauz_start*r_sre);  % Starting lag window index
        mz_end = floor(tauz_end*r_sre);  % Ending lag window index
        M = mz_end - mz_start + 1;
        % Equation 35 ECMA-418-2:2022
        phim_lztau = zeros(size(phim_lz_avg));
        phim_lztau(mz_start:mz_end, :) = phim_lz_avg(mz_start:mz_end, :) - mean(phim_lz_avg(mz_start:mz_end, :));  % lag-windowed, detrended ACF
        
        % Estimation of tonal loudness
        % ----------------------------
        % Section 6.2.5 Equation 36 ECMA-418-2:2022
        Phik_lztau = abs(fft(phim_lztau, 2*max(sz_b), 1));  % ACF spectrum in the lag window
        Phik_lztau(isnan(Phik_lztau)) = 0;
    
        % Section 6.2.5 Equation 37 ECMA-418-2:2022
        Nlz_tonal = phim_lz_avg(1, :);
        mask = 2*max(Phik_lztau, [], 1)/(M/2) <= phim_lz_avg(1, :);
        Nlz_tonal(mask) = 2*max(Phik_lztau(:, mask), [], 1)/(M/2);  % first estimation of specific loudness of tonal component in critical band
    
        % Section 6.2.5 Equations 38 & 39 ECMA-418-2:2022
        [~, kz_max] = max(Phik_lztau, [], 1);
        fz_ton = kz_max*(r_sre/(2*max(sz_b)));  % frequency of maximum tonal component in critical band

        % Section 6.2.7 Equation 41 ECMA-418-2:2022
        Nlz_signal = phim_lz_avg(1, :);  % specific loudness of complete band-pass signal in critical band
        
        % Resampling to common time basis Section 6.2.6 ECMA-418-2:2022
        if i_interp(zz) > 1
            % Note: use of interpolation function avoids rippling caused by
            % resample function, which disrupts specific loudness 
            % calculation for tonal and noise components
            l_n = size(phim_lz_avg, 2);
            x = linspace(1, l_n, l_n);
            xq = linspace(1, l_n, i_interp(zz)*l_n);
            Nlz_tonal = interp1(x, Nlz_tonal, xq);
            Nlz_signal = interp1(x, Nlz_signal, xq);
            fz_ton = interp1(x, fz_ton, xq);

        end

        % Remove end zero-padded samples Section 6.2.6 ECMA-418-2:2022
        % Note: in this part of the standard, the effect of the end
        % zero-padding to the input is removed from the processed signals.
        % There is no mention of realigning the processed signals to
        % compensate for the start zero-padding to the input. The
        % alternative terms below incorporate an amendment to the standard
        % to compensate for the start zero-padding, which results in
        % improved time alignment of the processed signals with the input.
        if ecma == true
            l_end = ceil(n_samples/r_sre*r_sd) + 1;  % Equation 40 ECMA-418-2:2022
            l_start = 1;  % start block
        else
            l_end = ceil((n_samples + n_zeross)/r_sre*r_sd) + 1;  % Equation 40 ECMA-418-2:2022 (edited to account for start zero-padding)
            l_start = floor(n_zeross/r_sre*r_sd) + 1;  % Additional term to remove start zero-padding lag
        end

        Nlz_tonal = Nlz_tonal(l_start:l_end);
        Nlz_signal = Nlz_signal(l_start:l_end);
        fz_ton = fz_ton(l_start:l_end);

        % Noise reduction Section 6.2.7 ECMA-418-2:2020
        % ---------------------------------------------
        SNRlz1 = Nlz_tonal./((Nlz_signal - Nlz_tonal) + 1e-12);  % Equation 42 ECMA-418-2:2022 signal-noise-ratio first approximation (ratio of tonal component loudness to non-tonal component loudness in critical band)
        Nlz_tonal = LowPass(Nlz_tonal, r_sd);  % Equation 43 ECMA-418-2:2022 low pass filtered specific loudness of non-tonal component in critical band
        SNRlz = LowPass(SNRlz1, r_sd);  % Equation 44 ECMA-418-2:2022 lowpass filtered SNR (improved estimation)
        gz = csz_b(zz)/(Fz(zz)^dsz_b(zz));  % Equation 46 ECMA-418-2:2022
        % Equation 45 ECMA-418-2:2022
        crit = exp(-alpha*((SNRlz/gz) - beta));
        nrlz = 1 - crit;  % sigmoidal weighting function
        nrlz(crit >= 1) = 0;
        Nlz_tonal = nrlz.*Nlz_tonal;  % Equation 47 ECMA-418-2:2022

        % Section 6.2.8 Equation 48 ECMA-418-2:2022
        Nlz_noise = LowPass(Nlz_signal, r_sd) - Nlz_tonal;  % specific loudness of non-tonal component in critical band
    
        % Store critical band results
        % ---------------------------
        SNRtz(zz, :, chan) = SNRlz1;  % specific time-dependent signal-noise-ratio in each critical band
        Ntz(zz, :, chan) = Nlz_signal;  % specific time-dependent loudness of signal in each critical band
        Ntz_ton(zz, :, chan) = Nlz_tonal;  % specific time-dependent loudness of tonal component in each critical band
        Ntz_noise(zz, :, chan) = Nlz_noise;   % specific time-dependent loudness of non-tonal component in each critical band
        ftz_ton(zz, :, chan) = fz_ton;  % time-dependent frequency of tonal component in each critical band
        
    end

    % Calculation of specific tonality
    % --------------------------------
    % Section 6.2.8 Equation 49 ECMA-418-2:2022
    SNRl = max(Ntz_ton, [], 1)./(1e-12 + sum(Ntz_noise, 1));  % loudness signal-noise-ratio
    
    % Section 6.2.8 Equation 50 ECMA-418-2:2022
    crit = exp(-A*(SNRl - B));
    ql = 1 - crit;  % sigmoidal scaling factor
    ql(crit >= 1) = 0;
    
    % Section 6.2.8 Equation 51 ECMA-418-2:2022#
    Ttz = c_T*ql.*Ntz_ton;  % time-dependent specific tonality
    
    % Calculation of time-averaged specific tonality Section 6.2.9 ECMA-418-2:2022
    for zz = 53:-1:1
        mask = Ttz(zz, :, chan) > 0.02;  % criterion Section 6.2.9 point 2
        mask(1:57) = 0;  % criterion Section 6.2.9 point 1
        % Section 6.2.9 Equation 53 ECMA-418-2:2022
        Tz(zz, 1, chan) = sum(Ttz(zz, mask, chan), 2)./(nnz(mask) + 1e-12);
        fz_ton(zz, 1, chan) = sum(ftz_ton(zz, mask, chan), 2)./(nnz(mask) + 1e-12);
    end

    % Calculation of total (non-specific) tonality Section 6.2.10
    % -----------------------------------------------------------
    % Further update can add the user input frequency range to determine
    % total tonality - not yet incorporated

    % Section 6.2.8 Equation 52 ECMA-418-2:2022
    t = (0:(size(Ttz, 2) - 1))/r_sd;  % time (s) corresponding with results output

    % Section 6.2.10 Equation 61 ECMA-418-2:2022
    [Tt(:, chan), zmax] = max(Ttz(:, :, chan), [], 1);  % Time-dependent total tonality
    for ll = size(ftz_ton, 2):-1:1
        ft_ton(1, ll, chan) = ftz_ton(zmax(ll), ll);
    end
    
    % Calculation of representative values Section 6.2.11 ECMA-418-2:2022
    % Time-averaged total tonality
    mask = Tt(:, chan) > 0.02;  % criterion Section 6.2.9 point 2
    mask(1:57) = 0;    % criterion Section 6.2.9 point 1
    % Section 6.2.11 Equation 63 ECMA-418-2:2022
    T(chan) = sum(Tt(mask, chan))/nnz(mask);  % Time-averaged total tonality (note: epsilon is not applied here, according to the standard)
    
    close(w)  % close waitbar

    % Plot figures
    % ------------
    if outplot
        % Plot results
        chan_lab = chans(chan);
        fig = figure;
        tile = tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, t, Fz, Ttz(:, :, chan), 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [t(1), t(end) + (t(2) - t(1))];
        ax1.YLim = [Fz(1), Fz(end)];
        ax1.CLim = [0, ceil(max(Tt(:, chan))*10)/10];
        ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
        ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]; 
        ax1.YScale = 'log';
        ax1.YLabel.String = 'Frequency, Hz';
        ax1.XLabel.String = 'Time, s';
        ax1.FontName =  'Arial';
        ax1.FontSize = 12;
        
        
        weightFilt = weightingFilter('A-weighting', r_sre);
        pA = weightFilt(p_re(:, chan));
        LA = 20*log10(rms(pA)/2e-5);
        axtitle = title(strcat(chan_lab, ' signal sound pressure level =', {' '}, num2str(round(LA,1)), "dB {\itL}_{Aeq}"),...
              'FontWeight', 'normal', 'FontName', 'Arial');
        colormap(turbo);
        h = colorbar;
        set(get(h,'label'),'string','Specific Tonality, Tu_{HMS}');
        
        ax2 = nexttile(2);
        plot(ax2, t, Tt(:, chan), 'r', 'LineWidth', 1);
        ax2.XLim = [t(1), t(end) + (t(2) - t(1))];
        ax2.YLim = [0, ceil(max(Tt(:, chan))*10)/10];
        ax2.XLabel.String = 'Time, s';
        ax2.YLabel.String = 'Tonality, Tu_{HMS}';
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.FontName =  'Arial';
        ax2.FontSize = 12;
    end

end

end

% LowPass Filter for Noise Reduction
function [y] = LowPass(x, r_s)

k = 3; % Footnote 21 ECMA-418-2:2022
e_i = [0, 1, 1]; % Footnote 21 ECMA-418-2:2022

% Footnote 20 ECMA-418-2:2022
tau = 1/32*6/7;

d = exp(-1/(r_s*tau)); % Section 5.1.4.2 ECMA-418-2:2022

% Feed-backward coefficients, Equation 14 ECMA-418-2:2022
m = 1:k;
a = [1, ((-d).^m).*arrayfun(@(m_) nchoosek(k, m_), m)];

% Feed-forward coefficients, Equation 15 ECMA-418-2:2022
m = 0:k-1;
i = 1:k-1;
b = (((1 - d)^k)./sum(e_i(i + 1).*(d.^i))).*(d.^m).*e_i;

% Recursive filter Equation 13 ECMA-418-2:2022
y = filter(b, a, x);

end

