function roughnessSHM = acousticSHMRoughness(p, sampleRateIn, axisN, soundField, waitBar, outPlot, binaural)
% roughnessSHM = acousticSHMRoughness(p, sampleRateIn, axisN, soundField, waitBar, outPlot, binaural)
%
% Returns roughness values according to ECMA-418-2:2025 (using the Sottek
% Hearing Model) for an input calibrated single mono or single stereo audio
% (sound pressure) time-series signal, p. For stereo signals, the binaural
% roughness can be calculated, and each channel is also analysed separately.
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
% soundField : keyword string (default: 'freeFrontal')
%              determines whether the 'freeFrontal' or 'diffuse' field stages
%              are applied in the outer-middle ear filter, or 'noOuter' uses
%              only the middle ear stage, or 'noEar' omits ear filtering.
%              Note: these last two options are beyond the scope of the
%              standard, but may be useful if recordings made using
%              artificial outer/middle ear are to be processed using the
%              specific recorded responses.
%
% waitBar : keyword string (default: true)
%           determines whether a progress bar displays during processing
%           (set waitBar to false for doing multi-file parallel calculations)
%
% outPlot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% binaural : Boolean true/false (default: true)
%            flag indicating whether to output combined binaural roughness
%            for stereo input signal.
%
% Returns
% -------
%
% roughnessSHM : structure
%                contains the output
%
% roughnessSHM contains the following outputs:
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
% If binaural=true, a corresponding set of outputs for the binaural
% roughness is also contained in roughnessSHM
%
% If outplot=true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall roughness, with the latter also indicating the time-aggregated
% value. A set of plots is returned for each input channel, with another
% set for the binaural roughness, if binaural=true. In that case, the
% indicated sound level corresponds with the channel with the highest sound
% level.
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
% Date created: 12/10/2023
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
        soundField (1, :) string {mustBeMember(soundField,...
                                               {'freeFrontal',...
                                                'diffuse', ...
                                                'noOuter', ...
                                                'noEar'})} = 'freeFrontal'
        waitBar {mustBeNumericOrLogical} = true
        outPlot {mustBeNumericOrLogical} = false
        binaural {mustBeNumericOrLogical} = true
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks
% Orient input matrix
if axisN == 2
    p = p.';
end

% Check the length of the input data (must be longer than 300 ms)
if size(p, 1) <=  300/1000*sampleRateIn
    error("Error: Input signal is too short along the specified axis to calculate roughness (must be longer than 300 ms)")
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

%% Define constants

signalT = size(p, 1)/sampleRateIn;  % duration of input signal
sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2025 [r_s]
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2025 [deltaf(f=0)]
c = 0.1618;  % Half-overlapping Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2025 [c]

dz = 0.5;  % critical band overlap [deltaz]
halfBark = 0.5:dz:26.5;  % half-overlapping critical band rate scale [z]
nBands = length(halfBark);  % number of bands
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025 [F(z)]

% Block and hop sizes Section 7.1.1 ECMA-418-2:2025
overlap = 0.75;  % block overlap proportion
blockSize = 16384;  % block size [s_b]
hopSize = (1 - overlap)*blockSize;  % hop size [s_h]

% Downsampled block and hop sizes Section 7.1.2 ECMA-418-2:2025
downSample = 32;  % downsampling factor
sampleRate1500 = sampleRate48k/downSample;
blockSize1500 = blockSize/downSample;
% hopSize1500 = (1 - overlap)*blockSize1500;
resDFT1500 = sampleRate1500/blockSize1500;  % DFT resolution (section 7.1.5.1) [deltaf]

% Modulation rate error correction values Table 8, Section 7.1.5.1
% ECMA-418-2:2025 [E(theta)]
errorCorrection = [0.0000, 0.0457, 0.0907, 0.1346, 0.1765, 0.2157, 0.2515,...
                   0.2828, 0.3084, 0.3269, 0.3364, 0.3348, 0.3188, 0.2844,...
                   0.2259, 0.1351, 0.0000];
errorCorrection = [errorCorrection, flip(-errorCorrection(1:end-1)), 0];

% High modulation rate roughness perceptual scaling function
% (section 7.1.5.2 ECMA-418-2:2025)
% Table 11 ECMA-418-2:2025 [r_1; r_2]
roughScaleParams = [0.3560, 0.8024;
                    0.8049, 0.9333];
roughScaleParams = [roughScaleParams(:, 1).*ones([2, sum(bandCentreFreqs < 1e3)]),...
                    roughScaleParams(:, 2).*ones([2, sum(bandCentreFreqs >= 1e3)])];
% Equation 84 ECMA-418-2:2025 [r_max(z)]
roughScale = 1./(1 + roughScaleParams(1, :).*abs(log2(bandCentreFreqs/1000)).^roughScaleParams(2, :));
roughScale = reshape(roughScale, [1, 1, nBands]);  % Note: this is to ease parallelised calculations

% High/low modulation rate roughness perceptual weighting function parameters
% (section 7.1.5.2 ECMA-418-2:2025)
% Equation 86 ECMA-418-2:2025 [f_max(z)]
modfreqMaxWeight = 72.6937*(1 - 1.1739*exp(-5.4583*bandCentreFreqs/1000));

% TODO for next update
%bandCentreFreqsWeight = max(ones(size(bandCentreFreqs)).*bandCentreFreqs(3), bandCentreFreqs);
%modfreqMaxWeight = 72.6937*(1 - 1.1739*exp(-5.4583*bandCentreFreqsWeight/1000));

% Equation 87 ECMA-418-2:2025 [q_1; q_2(z)]
roughHiWeightParams = [1.2822*ones(size(bandCentreFreqs));...
                       0.2471*ones(size(bandCentreFreqs))];
mask = bandCentreFreqs/1000 >= 2^-3.4253;
roughHiWeightParams(2, mask) = 0.2471 + 0.0129.*(log2(bandCentreFreqs(mask)/1000) + 3.4253).^2;
roughHiWeightParams = reshape(roughHiWeightParams, [2, 1, nBands]);  % Note: this is to ease parallelised calculations

% (section 7.1.5.4 ECMA-418-2:2025)
% Equation 96 ECMA-418-2:2025 [q_1; q_2(z)]
roughLoWeightParams = [0.7066*ones(size(bandCentreFreqs));...
                       1.0967 - 0.064.*log2(bandCentreFreqs/1000)];

% Output sample rate (section 7.1.7 ECMA-418-2:2025) [r_s50]
sampleRate50 = 50;

% Calibration constant
cal_R = 0.0180909;   % calibration factor in Section 7.1.7 Equation 104 ECMA-418-2:2025 [c_R]
%cal_Rx = 1/1.00123972659601;  % calibration adjustment factor
%cal_R*cal_Rx = 0.0180685; adjusted calibration value
cal_Rx = 1/1.0011565;  % calibration adjustment factor

% standardised epsilon
epsilon = 1e-12;

%% Signal processing

% Input pre-processing
% --------------------
if sampleRateIn ~= sampleRate48k  % Resample signal
    [p_re, ~] = shmResample(p, sampleRateIn);
else  % don't resample
    p_re = p;
end

% Input signal samples
n_samples = size(p_re, 1);

% Section 5.1.2 ECMA-418-2:2025 Fade in weighting and zero-padding
% (only the start is zero-padded)
pn = shmPreProc(p_re, blockSize, hopSize, true, false);

% Apply outer & middle ear filter
% -------------------------------
%
% Section 5.1.3.2 ECMA-418-2:2025 Outer and middle/inner ear signal filtering
pn_om = shmOutMidEarFilter(pn, soundField);

n_steps = 270;  % approximate number of calculation steps

% Loop through channels in file
% -----------------------------
for chan = chansIn:-1:1

    % Apply auditory filter bank
    % --------------------------
    
    if waitBar
        w = waitbar(0, "Initialising...");
        i_step = 1;
    
        waitbar(i_step/n_steps, w, 'Applying auditory filters...');
        i_step = i_step + 1;
    end % end of if branch for waitBar

    % Filter equalised signal using 53 1/2-overlapping Bark filters
    % according to Section 5.1.4.2 ECMA-418-2:2025
    pn_omz = shmAuditoryFiltBank(pn_om(:, chan), false);

    % Note: At this stage, typical computer RAM limits impose a need to loop
    % through the critical bands rather than continue with a parallelised
    % approach, until later downsampling is applied
    for zBand = nBands:-1:1
        % Segmentation into blocks
        % ------------------------
        if waitBar
            waitbar(i_step/n_steps, w, strcat("Calculating signal envelopes in 53 bands, ",...
                           num2str(zBand), " to go..."));...
            i_step = i_step + 1;
        end % end of if branch for waitBar

        % Section 5.1.5 ECMA-418-2:2025
        i_start = 1;
        [pn_lz, iBlocksOut] = shmSignalSegment(pn_omz(:, zBand), 1,...
                                               blockSize, overlap,...
                                               i_start, true);

        % Transformation into Loudness
        % ----------------------------
        if waitBar
            i_step = i_step + 1;
        end
        % Sections 5.1.6 to 5.1.9 ECMA-418-2:2025
        [~, bandBasisLoudness, ~] = shmBasisLoudness(pn_lz, bandCentreFreqs(zBand));
        basisLoudness(:, :, zBand) = bandBasisLoudness;
    
        % Envelope power spectral analysis
        % --------------------------------
        if waitBar
            i_step = i_step + 1;
        end
        % Sections 7.1.2 ECMA-418-2:2025
        % magnitude of Hilbert transform with downsample - Equation 65
        % [p(ntilde)_E,l,z]
        envelopes(:, :, zBand) = downsample(abs(hilbert(pn_lz)), downSample, 0);

    end  % end of for loop for obtaining low frequency signal envelopes

    % Note: With downsampled envelope signals, parallelised approach can continue

    % Section 7.1.3 equation 66 ECMA-418-2:2025 [Phi(k)_E,l,z]
    modSpectra = zeros(size(envelopes));
    envelopeWin = envelopes.*repmat(hann(blockSize1500, "periodic"), 1,...
                                    size(envelopes, 2), nBands)./sqrt(0.375);
    denom = max(basisLoudness, [], 3).*sum(envelopeWin.^2, 1);  % Equation 66 & 67
    mask = denom ~= 0;  % Equation 66 criteria for masking
    maskRep = repmat(mask, blockSize1500, 1 ,1);  % broadcast mask
    scaling = basisLoudness.^2./denom;  % Equation 66 factor
    scalingRep = repmat(scaling, blockSize1500, 1 ,1);  % broadcast scaling
    modSpectra(maskRep)...
        = reshape(scalingRep(maskRep),...
                  blockSize1500,...
                  [], nBands).*abs(fft(reshape(envelopeWin(maskRep),...
                                               blockSize1500, [], nBands))).^2;

    % Envelope noise reduction
    % ------------------------
    % section 7.1.4 ECMA-418-2:2025
    modSpectraAvg = modSpectra;
    modSpectraAvg(:, :, 2:end-1) = movmean(modSpectra, [1, 1], 3, 'Endpoints', 'discard');
    
    modSpectraAvgSum = sum(modSpectraAvg, 3);  % Equation 68 [s(l,k)]

    % Equation 71 ECMA-418-2:2025 [wtilde(l,k)]
    clipWeight = 0.0856.*modSpectraAvgSum(1:size(modSpectraAvg, 1)/2 + 1, :)...
                 ./(median(modSpectraAvgSum(3:size(modSpectraAvg, 1)/2, :), 1) + 1e-10)...
                 .*transpose(min(max(0.1891.*exp(0.012.*(0:size(modSpectraAvg, 1)/2)), 0), 1));

    % Equation 70 ECMA-418-2:2025 [w(l,k)]
    weightingFactor1 = zeros(size(modSpectraAvgSum(1:257, :, :)));
    mask = clipWeight >= 0.05*max(clipWeight(3:256, :), [], 1);
    weightingFactor1(mask) = min(max(clipWeight(mask) - 0.1407, 0), 1);
    weightingFactor = [weightingFactor1;
                       flipud(weightingFactor1(2:256, :))];

    % Calculate noise-reduced, scaled, weighted modulation power spectra
    modWeightSpectraAvg = modSpectraAvg.*weightingFactor; % Equation 69 [Phihat(k)_E,l,z]

    % Spectral weighting
    % ------------------
    % Section 7.1.5 ECMA-418-2:2025
    % theta used in equation 79, including additional index for
    % errorCorrection terms from table 10
    theta = 0:1:33;
    mlabIndex = 1;  % term used to compensate for MATLAB 1-indexing
    nBlocks = size(modWeightSpectraAvg, 2);
    modAmp = zeros(10, nBlocks, nBands);
    modRate = zeros(10, nBlocks, nBands);

    for zBand = nBands:-1:1
        if waitBar
            waitbar(i_step/n_steps, w, strcat("Calculating spectral weightings in 53 bands, ",...
                           num2str(zBand), " to go..."));...
            i_step = i_step + 1;
        end % end of if branch for waitBar

        % Section 7.1.5.1 ECMA-418-2:2025
        for lBlock = nBlocks:-1:1
            % identify peaks in each block (for each band)
            startIdx = 2;
            endIdx = 255;
            [PhiPks, kLocs, ~, proms] = findpeaks(modWeightSpectraAvg(startIdx + mlabIndex:endIdx + mlabIndex,...
                                                                      lBlock,...
                                                                      zBand));

            % reindex kLocs to match spectral start index used in findpeaks
            % for indexing into modulation spectra matrices
            kLocs = kLocs + startIdx;

            % we can only have peaks at k = 3:254
            mask = ismember(kLocs, 3:254);
            kLocs = kLocs(mask);
            PhiPks = PhiPks(mask);

            % consider 10 highest prominence peaks only
            if length(proms) > 10
                [promsSorted, iiSort] = sort(proms, 'descend');
                mask = proms >= promsSorted(10);
                
                % if branch to deal with duplicated peak prominences
                if sum(mask) > 10
                   mask = mask(iiSort <= 10); 
                end  % end of if branch for duplicated peak prominences

                PhiPks = PhiPks(mask);
                kLocs = kLocs(mask);

            end  % end of if branch to select 10 highest prominence peaks

            % consider peaks meeting criterion
            if ~isempty(PhiPks)
                mask = PhiPks > 0.05*max(PhiPks);  % Equation 72 criterion
                PhiPks = PhiPks(mask);  % [Phihat(k_p,i(l,z))]
                kLocs = kLocs(mask);
                % loop over peaks to obtain modulation rates
                for iPeak = length(PhiPks):-1:1
                    % Equation 74 ECMA-418-2:2025
                    % Note: here, the kLoc values are used as indices for
                    % the modulation spectral matrix, so MATLAB indexing
                    % is correctly addressed (see Equation 75 below)
                    % [Phihat_E,l,z]
                    modAmpMat = [modWeightSpectraAvg(kLocs(iPeak) - 1, lBlock, zBand);
                                 modWeightSpectraAvg(kLocs(iPeak), lBlock, zBand);
                                 modWeightSpectraAvg(kLocs(iPeak) + 1, lBlock, zBand)];
                    
                    % Equation 82 [A_i(l,z)]
                    modAmp(iPeak, lBlock, zBand) = sum(modAmpMat);

                    % Equation 75 ECMA-418-2:2025
                    % Note: because the kLoc index values are used directly
                    % in the calculation, MATLAB indexing needs to be
                    % compensated for by subtracting 1 from kLocs
                    % [K]
                    modIndexMat = [(kLocs(iPeak) - mlabIndex - 1)^2, kLocs(iPeak) - mlabIndex - 1, 1;
                                   (kLocs(iPeak) - mlabIndex)^2, kLocs(iPeak) - mlabIndex, 1;
                                   (kLocs(iPeak) - mlabIndex + 1)^2, kLocs(iPeak) - mlabIndex + 1, 1];

                    coeffVec = modIndexMat\modAmpMat;  % Equation 73 solution [C]

                    % Equation 76 ECMA-418-2:2025 [ftilde_p,i(l,z)]
                    modRateEst = -(coeffVec(2)/(2*coeffVec(1)))*resDFT1500;

                    % Equation 79 ECMA-418-2:2025 [beta(theta)]
                    errorBeta = (floor(modRateEst/resDFT1500) + theta(1:33)/32)*resDFT1500...
                                - (modRateEst + errorCorrection(theta(1:33) + mlabIndex)); % compensated theta value for MATLAB-indexing

                    % Equation 80 ECMA-418-2:2025 [theta_min]
                    [~, i_minError] = min(abs(errorBeta));
                    thetaMinError = theta(i_minError);  % the result here is 0-indexed

                    % Equation 81 ECMA-418-2:2025 [theta_corr]
                    if thetaMinError > 0 && errorBeta(i_minError)*errorBeta(i_minError - 1) < 0 % (0-indexed)
                        thetaCorr = thetaMinError;  % 0-indexed
                    else
                        thetaCorr = thetaMinError + 1;  % 0-indexed
                    end  % end of eq 81 if-branch

                    % Equation 78 ECMA-418-2:2025
                    % thetaCorr is 0-indexed so needs adjusting when
                    % indexing
                    % [rho(ftilde_p,i(l,z))]
                    biasAdjust = errorCorrection(thetaCorr + mlabIndex - 1)...
                                 - (errorCorrection(thetaCorr + mlabIndex)...
                                    - errorCorrection(thetaCorr + mlabIndex - 1))...
                                    *errorBeta(thetaCorr + mlabIndex - 1)...
                                    /(errorBeta(thetaCorr + mlabIndex)...
                                      - errorBeta(thetaCorr + mlabIndex - 1));

                    % Equation 77 ECMA-418-2:2025 [f_p,i(l,z)]
                    modRate(iPeak, lBlock, zBand) = modRateEst + biasAdjust;

                end  % end of for loop over peaks in block per band
            end  % end of if branch for detected peaks in modulation spectrum
        end  % end of for loop over blocks for peak detection      
    end  % end of for loop over bands for modulation spectral weighting

    % Section 7.1.5.2 ECMA-418-2:2025 - Weighting for high modulation rates
    % Equation 85 [G_l,z,i(f_p,i(l,z))]
    roughHiWeight = shmRoughWeight(modRate,...
                                   reshape(modfreqMaxWeight, [1, 1, nBands]),...
                                   roughHiWeightParams);

    % Equation 83 [Atilde_i(l,z)]
    modAmpHiWeight = modAmp.*roughScale;
    mask = modRate <= resDFT1500;
    modAmpHiWeight(mask) = 0;
    mask = modRate > permute(repmat(modfreqMaxWeight, 1, 1, 10), [3, 1, 2]);
    modAmpHiWeight(mask) = modAmpHiWeight(mask).*roughHiWeight(mask);


    % Section 7.1.5.3 ECMA-418-2:2025 - Estimation of fundamental modulation rate
    % TODO: replace the loop approach with a parallelised approach!
    % matrix initialisation to ensure zero rates do not cause missing bands in output
    modFundRate = zeros([nBlocks, nBands]);
    modMaxWeight = zeros([10, nBlocks, nBands]);
    for zBand = nBands:-1:1
        if waitBar
            waitbar(i_step/n_steps, w, strcat("Calculating modulation rates in 53 bands, ",...
                    num2str(zBand), " to go..."));...
            i_step = i_step + 1;
        end % end of if branch for waitBar

        for lBlock = nBlocks:-1:1
            % Proceed with rate detection if non-zero modulation rates
            if max(modRate(:, lBlock, zBand)) > 0
                modRateForLoop = modRate(modRate(:, lBlock, zBand) > 0,...
                                         lBlock, zBand);

                nPeaks = length(modRateForLoop);

                % initialise empty cell array for equation 90
                indSetiPeak = {};
                % initialise empty matrix for equation 91
                harmCompEnergy = double.empty(nPeaks, 0);

                for iPeak = nPeaks:-1:1
                    % Equation 88 [R_i_0(i)]
                    modRateRatio = round(modRateForLoop/modRateForLoop(iPeak));
                    [uniqRatios, startGroupInds, uniqGroupInds] = unique(modRateRatio);
                    countDupes = accumarray(uniqGroupInds, 1);

                    % add any non-duplicated ratio indices
                    testIndices = zeros([10, 1]);
                    if ~isempty(startGroupInds(countDupes==1))
                        testIndices(1:length(startGroupInds(countDupes==1))) = startGroupInds(countDupes==1);
                    end

                    % loop over duplicated values to select single
                    % index
                    if max(countDupes) > 1
                        dupeRatioVals = uniqRatios(countDupes > 1);
                        for jDupe = length(dupeRatioVals):-1:1

                            % Equation 89 [i]
                            dupeGroupInds = find(modRateRatio == dupeRatioVals(jDupe));
                            testDupe = abs(modRateForLoop(dupeGroupInds)...
                                           ./(modRateRatio(dupeGroupInds)...
                                            *modRateForLoop(iPeak)) - 1);

                            % discard if all inf
                            if ~all(isinf(testDupe))
                                [~, testDupeMin] = min(testDupe);
                                % append selected index
                                testIndices(length(startGroupInds(countDupes==1)) + jDupe) = dupeGroupInds(testDupeMin);
                            end  % end of if branch for all inf
                        end  % end of for loop over duplicated ratios
                    end  % end of if branch for duplicated ratios

                    % discard zero indices
                    testIndices = testIndices(testIndices > 0);

                    % Equation 90 [I_i_0]
                    harmComplexTest = abs(modRateForLoop(testIndices)./(modRateRatio(testIndices)*modRateForLoop(iPeak)) - 1);
                    indSetiPeak{iPeak} = testIndices(harmComplexTest < 0.04);

                    % Equation 91 [E_i_0]
                    harmCompEnergy(iPeak) = sum(modAmpHiWeight(indSetiPeak{iPeak}, lBlock, zBand));

                end

                [~, iMaxEnergy] = max(harmCompEnergy);
                indSetMax = indSetiPeak{iMaxEnergy};
                modFundRate(lBlock, zBand) = modRateForLoop(iMaxEnergy);
                % Equation 94 [i_peak]
                [~, iPeakAmp] = max(modAmpHiWeight(indSetMax, lBlock, zBand));
                iPeak = indSetMax(iPeakAmp);
    
                % Equation 93 [w_peak]
                gravityWeight = 1 + 0.1*abs(sum(modRateForLoop(indSetMax)...
                                            .*modAmpHiWeight(indSetMax, lBlock, zBand))...
                                            /sum(modAmpHiWeight(indSetMax, lBlock, zBand) + epsilon)...
                                            - modRateForLoop(iPeak)).^0.749;

                % Equation 92 [Ahat(i)]
                modMaxWeight(indSetMax,...
                             lBlock,...
                             zBand) = gravityWeight.*modAmpHiWeight(indSetMax,...
                                                                    lBlock,...
                                                                    zBand);

            end  % end of if branch for non-zero modulation rates
        end  % end of for loop over blocks
    end  % end of for loop over bands

    % Equation 95 [A(l,z)]
    roughLoWeight = shmRoughWeight(modFundRate, modfreqMaxWeight, roughLoWeightParams);
    modMaxWeightSum = squeeze(sum(modMaxWeight, 1));
    modMaxLoWeight = squeeze(sum(permute(repmat(roughLoWeight, 1, 1, 10), [3, 1, 2]).*modMaxWeight, 1));
    mask = modFundRate <= resDFT1500;
    modMaxLoWeight(mask) = 0;
    mask = modFundRate > modfreqMaxWeight;
    modMaxLoWeight(mask) = modMaxWeightSum(mask);
    modAmpMax = modMaxLoWeight;
    modAmpMax(modAmpMax < 0.074376) = 0;

    % Time-dependent specific roughness
    % ---------------------------------
    % Section 7.1.7 ECMA-418-2:2025

    % interpolation to 50 Hz sampling rate
    % Section 7.1.7 Equation 103 [l_50,end]
    l_50Last = floor(n_samples/sampleRate48k*sampleRate50) + 1;
    x = (iBlocksOut - 1)/sampleRate48k;
    xq = linspace(0, signalT, l_50Last);
    for zBand = nBands:-1:1
        specRoughEst(:, zBand) = pchip(x, modAmpMax(:, zBand), xq);
    end  % end of for loop for interpolation
    specRoughEst(specRoughEst < 0) = 0;  % [R'_est(l_50,z)]

    % Section 7.1.7 Equation 107 [Rtilde'_est(l_50)]
    specRoughEstRMS = rms(specRoughEst, 2);

    % Section 7.1.7 Equation 108 [Rbar'_est(l_50)]
    specRoughEstAvg = mean(specRoughEst, 2);

    % Section 7.1.7 Equation 106 [B(l_50)]
    Bl50 = zeros(size(specRoughEstAvg));
    mask = specRoughEstAvg ~= 0;
    Bl50(mask) = specRoughEstRMS(mask)./specRoughEstAvg(mask);

    % Section 7.1.7 Equation 105 [E(l_50)]
    El50 = (0.95555 - 0.58449)*(tanh(1.6407*(Bl50 - 2.5804)) + 1)*0.5 + 0.58449;

    % Section 7.1.7 Equation 104 [Rhat'(l_50,z)]
    specRoughEstTform = cal_R*cal_Rx*(specRoughEst.^El50);

    % Section 7.1.7 Equation 109-110 [R'(l_50,z)]
    riseTime = 0.0625;
    fallTime = 0.5;
    specRoughness(:, :, chan) = shmRoughLowPass(specRoughEstTform, sampleRate50, ...
                                                riseTime, fallTime);

    if waitBar
        close(w)  % close waitbar
    end

end  % end of for loop over channels

% Binaural roughness
% Section 7.1.11 ECMA-418-2:2025 [R'_B(l_50,z)]
if chansIn == 2 && binaural
    % Equation 112
    specRoughness(:, :, 3) = sqrt(sum(specRoughness(:, :, 1:2).^2, 3)/2);
    chansOut = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Binaural"];
else
    chansOut = chansIn;  % assign number of output channels
end

% Section 7.1.8 ECMA-418-2:2025
% Time-averaged specific roughness [R'(z)]
specRoughnessAvg = mean(specRoughness(17:end, :, :), 1);

% Section 7.1.9 ECMA-418-2:2025
% Time-dependent roughness Equation 111 [R(l_50)]
% Discard singleton dimensions
if chansOut == 1
    roughnessTDep = sum(specRoughness.*dz, 2);
    specRoughnessAvg = transpose(specRoughnessAvg);
else
    roughnessTDep = squeeze(sum(specRoughness.*dz, 2));
    specRoughnessAvg = squeeze(specRoughnessAvg);
end

% Section 7.1.10 ECMA-418-2:2025
% Overall roughness [R]
roughness90Pc = prctile(roughnessTDep(17:end, :), 90, 1);

% time (s) corresponding with results output [t]
timeOut = (0:(size(specRoughness, 1) - 1))/sampleRate50;

%% Output plotting

if outPlot
    % Plot figures
    % ------------
    for chan = chansOut:-1:1
        % Plot results
        fig = figure;
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, timeOut, bandCentreFreqs, permute(specRoughness(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        ax1.CLim = [0, ceil(max(specRoughness(:, :, chan), [], 'all')*500)/500];
        ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
        ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                          "8k", "16k"];
        ax1.YScale = 'log';
        ax1.YLabel.String = 'Frequency, Hz';
        ax1.XLabel.String = 'Time, s';
        ax1.FontName =  'Arial';
        ax1.FontSize = 12;
        cmap_inferno = load('cmap_inferno.txt');
        colormap(cmap_inferno);
        h = colorbar;
        set(get(h,'label'),'string', {'Specific roughness,'; 'asper_{SHM}/Bark_{SHM}'});        
        chan_lab = chans(chan);

        % Create A-weighting filter
        weightFilt = weightingFilter('A-weighting', sampleRateIn);
        % Filter signal to determine A-weighted time-averaged level
        if chan == 3
            pA = weightFilt(p);
            LAeq2 = 20*log10(rms(pA, 1)/2e-5);
            % take the higher channel level as representative (PD ISO/TS
            % 12913-3:2019 Annex D)
            [LAeq, LR] = max(LAeq2);
            % if branch to identify which channel is higher
            if LR == 1
                whichEar = ' left ear';
            else
                whichEar = ' right ear';
            end  % end of if branch

            chan_lab = chan_lab + whichEar;

        else
            pA = weightFilt(p(:, chan));
            LAeq = 20*log10(rms(pA)/2e-5);
        end
        
        title(strcat(chan_lab,...
                     ' signal sound pressure level =', {' '},...
                     num2str(round(LAeq,1)), "dB {\itL}_{Aeq}"),...
                     'FontWeight', 'normal', 'FontName', 'Arial');

        ax2 = nexttile(2);
        plot(ax2, timeOut, roughness90Pc(1, chan)*ones(size(timeOut)), 'color',...
             cmap_inferno(34, :), 'LineWidth', 1, 'DisplayName', "90th" + string(newline) + "percentile");
        hold on
        plot(ax2, timeOut, roughnessTDep(:, chan), 'color', cmap_inferno(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        if max(roughnessTDep(:, chan)) > 0
            ax2.YLim = [0, 1.1*ceil(max(roughnessTDep(:, chan))*10)/10];
        end
        ax2.XLabel.String = 'Time, s';
        ax2.YLabel.String = 'Roughness, asper_{SHM}';
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.GridAlpha = 0.075;
        ax2.GridLineStyle = '--';
        ax2.GridLineWidth = 0.25;
        ax2.FontName = 'Arial';
        ax2.FontSize = 12;
        lgd = legend('Location', 'eastoutside', 'FontSize', 8);
        lgd.Title.String = "Overall";
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting if outplot true

%% Output assignment

% Assign outputs to structure
if chansOut == 3
    roughnessSHM.specRoughness = specRoughness(:, :, 1:2);
    roughnessSHM.specRoughnessAvg = specRoughnessAvg(:, 1:2);
    roughnessSHM.roughnessTDep = roughnessTDep(:, 1:2);
    roughnessSHM.roughness90Pc = roughness90Pc(:, 1:2);
    roughnessSHM.specRoughnessBin = specRoughness(:, :, 3);
    roughnessSHM.specRoughnessAvgBin = specRoughnessAvg(:, 3);
    roughnessSHM.roughnessTDepBin = roughnessTDep(:, 3);
    roughnessSHM.roughness90PcBin = roughness90Pc(:, 3);
    roughnessSHM.bandCentreFreqs = bandCentreFreqs;
    roughnessSHM.timeOut = timeOut;
    roughnessSHM.soundField = soundField;
else
    roughnessSHM.specRoughness = specRoughness;
    roughnessSHM.specRoughnessAvg = specRoughnessAvg;
    roughnessSHM.roughnessTDep = roughnessTDep;
    roughnessSHM.roughness90Pc = roughness90Pc;
    roughnessSHM.bandCentreFreqs = bandCentreFreqs;
    roughnessSHM.timeOut = timeOut;
    roughnessSHM.soundField = soundField;
end

% end of function
