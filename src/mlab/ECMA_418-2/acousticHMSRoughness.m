function roughnessHMS = acousticHMSRoughness(p, sampleRatein, axisn, outplot, binaural)
% roughnessHMS = acousticHMSRoughness(p, sampleRatein, axisn, outplot, binaural)
%
% Returns roughness values **and frequencies** according to ECMA-418-2:2022
% (using the Hearing Model of Sottek) for an input calibrated single mono
% or single stereo audio (sound pressure) time-series signal, p.
%
% Inputs
% ------
% p : vector or 2D matrix
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% sampleRatein : integer
%                the sample rate (frequency) of the input signal(s)
%
% axisn : integer (1 or 2, default: 1)
%         the time axis along which to calculate the tonality
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% binaural : Boolean true/false (default: true)
%            flag indicating whether to output binaural roughness for stereo
%            input signal.
% 
% Returns
% -------
%
% roughnessHMS : structure
%                contains the output
%
% roughnessHMS contains the following outputs:
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
% If binaural=true, a corresponding set of outputs for the binaural
% roughness is also contained in roughnessHMS
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
% Date last modified: 13/08/2024
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
        sampleRatein (1, 1) double {mustBePositive, mustBeInteger}
        axisn (1, 1) {mustBeInteger, mustBeInRange(axisn, 1, 2)} = 1
        outplot {mustBeNumericOrLogical} = false
        binaural {mustBeNumericOrLogical} = true
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))

%% Input checks
% Orient input matrix
if axisn == 2
    p = p.';
end

% Check the length of the input data (must be longer than 300 ms)
if size(p, 1) <=  300/1000*sampleRatein
    error("Error: Input signal is too short to calculate roughness (must be longer than 300 ms)")
end

% Check the channel number of the input data
if size(p, 2) > 2
    error("Error: Input signal comprises more than two channels")
else
    inchans = size(p, 2);
    if inchans == 2
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

%% Define constants

signalT = size(p, 1)/sampleRatein;  % duration of input signal
sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2022
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2022
c = 0.1618;  % Half-Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2022

dz = 0.5;  % critical band resolution
halfBark = 0.5:dz:26.5;  % half-critical band rate scale
nBands = length(halfBark);  % number of bands
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2022
dfz = sqrt(deltaFreq0^2 + (c*bandCentreFreqs).^2);  % Section 5.1.4.1 Equation 10 ECMA-418-2:2022

% Block and hop sizes Section 7.1.1 ECMA-418-2:2022
overlap = 0.75;  % block overlap proportion
blockSize = 16384;  % block size
hopSize = (1 - overlap)*blockSize;  % hop size

% Downsampled block and hop sizes Section 7.1.2 ECMA-418-2:2022
downSample = 32;  % downsampling factor
sampleRate1500 = sampleRate48k/downSample;
blockSize1500 = blockSize/downSample;
hopSize1500 = (1 - overlap)*blockSize1500;
resDFT1500 = sampleRate1500/blockSize1500;  % DFT resolution (section 7.1.5.1)

% Modulation rate error correction values Table 8, Section 7.1.5.1 ECMA-418-2:2022
errorCorrection = [0.0000, 0.0457, 0.0907, 0.1346, 0.1765, 0.2157, 0.2515,...
                   0.2828, 0.3084, 0.3269, 0.3364, 0.3348, 0.3188, 0.2844,...
                   0.2259, 0.1351, 0.0000];
errorCorrection = [errorCorrection, flip(-errorCorrection(1:end-1)), 0];

% High modulation rate roughness perceptual scaling function
% (section 7.1.5.2 ECMA-418-2:2022)
% Table 11 ECMA-418-2:2022
roughScaleParams = [0.3560, 0.8024;
                    0.8049, 0.9333];
roughScaleParams = [roughScaleParams(:, 1).*ones([2, sum(bandCentreFreqs < 1e3)]),...
                    roughScaleParams(:, 2).*ones([2, sum(bandCentreFreqs >= 1e3)])];
% Equation 84 ECMA-418-2:2022
roughScale = 1./(1 + roughScaleParams(1, :).*abs(log2(bandCentreFreqs/1000)).^roughScaleParams(2, :));
roughScale = reshape(roughScale, [1, 1, nBands]);  % Note: this is to ease parallelised calculations

% High modulation rate roughness perceptual weighting function parameters
% (section 7.1.5.2 ECMA-418-2:2022)
% Equation 86 ECMA-418-2:2022
modfreqMaxWeight = 72.6937*(1 - 1.1739*exp(-5.4583*bandCentreFreqs/1000));

% Equation 87 ECMA-418-2:2022
roughHiWeightParams = [1.2822*ones(size(bandCentreFreqs));...
                       0.2471*ones(size(bandCentreFreqs))];
mask = bandCentreFreqs/1000 >= 2^-3.4253;
roughHiWeightParams(2, mask) = 0.2471 + 0.0129.*(log2(bandCentreFreqs(mask)/1000) + 3.4253).^2;
roughHiWeightParams = reshape(roughHiWeightParams, [2, 1, nBands]);  % Note: this is to ease parallelised calculations

% Equation 96 ECMA-418-2:2022
roughLoWeightParams = [0.7066*ones(size(bandCentreFreqs));...
                       1.0967 - 0.064.*log2(bandCentreFreqs/1000)];

% Output sample rate (section 7.1.7 ECMA-418-2:2022)
sampleRate50 = 50;

% Calibration constant
c_R = 0.0180685;
c_Rx = 1;  % calibration adjustment factor

%% Signal processing

% Input pre-processing
% --------------------
if sampleRatein ~= sampleRate48k  % Resample signal
    [p_re, ~] = hmSResample(p, sampleRatein);
else  % don't resample
    p_re = p;
end

% Input signal samples
n_samples = size(p_re, 1);

% Section 5.1.2 ECMA-418-2:2022 Fade in weighting and zero-padding
% (only the start is zero-padded)
pn = hmSPreProc(p_re, max(blockSize), max(hopSize), true, false);

% Apply outer & middle ear filter
% -------------------------------
%
% Section 5.1.3.2 ECMA-418-2:2022 Outer and middle/inner ear signal filtering
pn_om = hmSOutMidEarFilter(pn);

n_steps = 250;  % approximate number of calculation steps

% Loop through channels in file
% -----------------------------
for chan = size(pn_om, 2):-1:1
    w = waitbar(0, "Initialising...");
    i_step = 1;

    % Apply auditory filter bank
    % --------------------------
    waitbar(i_step/n_steps, w, 'Applying auditory filters...');
    i_step = i_step + 1;

    % Filter equalised signal using 53 1/2Bark ERB filters according to 
    % Section 5.1.4.2 ECMA-418-2:2022
    pn_omz = hmSAuditoryFiltBank(pn_om(:, chan), false);

    % Note: At this stage, typical computer RAM limits impose a need to loop
    % through the critical bands rather than continue with a parallelised
    % approach, until later downsampling is applied
    for zBand = nBands:-1:1
        % Segmentation into blocks
        % ------------------------
        waitbar(i_step/n_steps, w, "Calculating signal envelopes in 53 bands, ",...
                       num2str(zBand), " to go...");...
        i_step = i_step + 1;

        % Section 5.1.5 ECMA-418-2:2022
        i_start = 1;
        [pn_lz, iBlocks] = hmSsignalSegment(pn_omz(:, zBand), 1, blockSize, overlap,...
                                            i_start, true);

        % Transformation into Loudness
        % ----------------------------
        i_step = i_step + 1;
        % Sections 5.1.6 to 5.1.9 ECMA-418-2:2022
        [~, bandBasisLoudness, ~] = hmSBasisLoudness(pn_lz, bandCentreFreqs(zBand));
        basisLoudness(:, :, zBand) = bandBasisLoudness;
    
        % Envelope power spectral analysis
        % --------------------------------
        i_step = i_step + 1;
        % Sections 7.1.2 ECMA-418-2:2022
        % magnitude of Hilbert transform with downsample - Equation 65
        % [p(ntilde)_E,l,z]
        % Notefigure; imagesc: prefiltering is not needed because the output from the
        % Hilbert transform is a form of low-pass filtering
        envelopes(:, :, zBand) = downsample(abs(hilbert(pn_lz)), downSample, 0);

    end  % end of for loop for obtaining low frequency signal envelopes

    % Note: With downsampled envelope signals, parallelised approach can continue

    % Section 7.1.3 equation 66 ECMA-418-2:2022 [Phi(k)_E,l,z]
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
    % section 7.1.4 ECMA-418-2:2022
    modSpectraAvg = movmean(modSpectra, [1, 1], 3, 'Endpoints', 'shrink');
    modSpectraAvgSum = sum(modSpectraAvg, 3);  % Equation 68 [s(l,k)]

    % Equation 71 ECMA-418-2:2022 [wtilde(l,k)]
    clipWeight = 0.0856.*modSpectraAvgSum(1:size(modSpectraAvg, 1)/2, :)...
                 ./(median(modSpectraAvgSum(3:size(modSpectraAvg, 1)/2, :), 1) + 1e-10)...
                 .*transpose(min(max(0.1891.*exp(0.012.*(0:size(modSpectraAvg, 1)/2 - 1)), 0), 1));

    % Equation 70 ECMA-418-2:2022
    weightingFactor1 = zeros(size(modSpectraAvgSum(1:256, :, :)));
    mask = clipWeight >= 0.05*max(clipWeight(3:256, :), [], 1);
    weightingFactor1(mask) = min(max(clipWeight(mask) - 0.1407, 0), 1);
    weightingFactor = [weightingFactor1;
                       flipud(weightingFactor1)];

    % Calculate noise-reduced, scaled, weighted modulation power spectra
    modWeightSpectraAvg = modSpectraAvg.*weightingFactor; % Equation 69 [Phihat(k)_E,l,z]

    % Spectral weighting
    % ------------------
    % Section 7.1.5 ECMA-418-2:2022
    n_blocks = size(modWeightSpectraAvg, 2);
    for zBand = nBands:-1:1
        waitbar(i_step/n_steps, w, "Calculating spectral weightings in 53 bands, ",...
                       num2str(zBand), " to go...");...
        i_step = i_step + 1;

        % Section 7.1.5.1 ECMA-418-2:2022
        for ll = n_blocks:-1:1
            % assign zero arrays for ten greatest spectral maxima amplitudes
            % and modulation rates in each block
            modAmpBlock = zeros(10, 1);
            modRateBlock = zeros(10, 1);
            % identify peaks in each block (for each band)
            [PhiPks, kLocs, ~, proms] = findpeaks(modWeightSpectraAvg(3:256,...
                                                                      ll,...
                                                                      zBand));

            % reindex locs to match spectral start index used in findpeaks
            kLocs = kLocs + 2;

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
                PhiPks = PhiPks(mask);
                kLocs = kLocs(mask);
                % loop over peaks to obtain modulation rates
                for ii = length(PhiPks):-1:1
                    % Equation 74 ECMA-418-2:2022
                    modAmpMat = [modWeightSpectraAvg(kLocs(ii) - 1, ll, zBand);
                               modWeightSpectraAvg(kLocs(ii), ll, zBand);
                               modWeightSpectraAvg(kLocs(ii) + 1, ll, zBand)];
                    
                    modAmpBlock(ii) = sum(modAmpMat);

                    % Equation 75 ECMA-418-2:2022
                    % Note: because the kLoc values are used directly in
                    % the calculation, MATLAB indexing needs to be
                    % accounted for
                    modIndexMat = [(kLocs(ii) - 2)^2, kLocs(ii) - 2, 1;
                                   (kLocs(ii) - 1)^2, kLocs(ii) - 1,     1;
                                   (kLocs(ii))^2, kLocs(ii), 1];

                    coeffVec = modIndexMat\modAmpMat;  % Equation 73 solution

                    % Equation 76 ECMA-418-2:2022
                    modRateEst = -(coeffVec(2)/(2*coeffVec(1)))*resDFT1500;

                    % Equation 79 ECMA-418-2:2022
                    theta = (0:1:33) + 1;
                    errorBeta = (floor(modRateEst/resDFT1500) + (theta - 1)/32)*resDFT1500...
                                - (modRateEst + errorCorrection);

                    % Equation 80 ECMA-418-2:2022
                    [~, i_minError] = min(abs(errorBeta));
                    thetaMinError = theta(i_minError);

                    % Equation 81 ECMA-418-2:2022
                    if  thetaMinError > 1 && (floor(modRateEst/resDFT1500) + (thetaMinError - 1)/32)*resDFT1500...
                                - (modRateEst + errorCorrection(thetaMinError)) < 0
                        thetaCorr = thetaMinError;
                    else
                        thetaCorr = thetaMinError + 1;
                    end  % end of eq 81 if branch

                    % Equation 78 ECMA-418-2:2022
                    biasAdjust = errorCorrection(thetaCorr)...
                                 - (errorCorrection(thetaCorr - 1)...
                                    - errorCorrection(thetaCorr - 1))...
                                    *(errorBeta(thetaCorr - 1)/(errorBeta(thetaCorr) - errorBeta(thetaCorr - 1)));

                    % Equation 77 ECMA-418-2:2022
                    modRateBlock(ii) = modRateEst + biasAdjust;

                end  % end of for loop over peaks in block per band
            end  % end of if branch for detected peaks in modulation spectrum

            % collect modulation peak amplitudes and frequencies for all
            % blocks in each band
            modAmpBand(:, ll) = modAmpBlock;
            modRateBand(:, ll) = modRateBlock;
        end  % end of for loop over blocks for peak detection
        
        % collect modulation peak amplitudes and frequencies for all bands
        modAmp(:, :, zBand) = modAmpBand;
        modRate(:, :, zBand) = modRateBand;
        
    end  % end of for loop over bands for modulation spectral weighting

    % Section 7.1.5.2 ECMA-418-2:2022 - Weighting for high modulation rates
    % Equation 85
    roughHiWeight = hmSRoughWeight(modRate,...
                                   reshape(modfreqMaxWeight, [1, 1, nBands]),...
                                   roughHiWeightParams);

    % Equation 83
    modHiWeight = modAmp.*roughScale;
    mask = modRate >= permute(repmat(modfreqMaxWeight, 1, 1, 10), [3, 1, 2]);
    modHiWeight(mask) = modHiWeight(mask).*roughHiWeight(mask);

    % Section 7.1.5.3 ECMA-418-2:2022 - Estimation of fundamental modulation rate

    % the loop approach - this section is based on a translation of
    % _estimate_fund_mod_rate.py from the Python MoSQITo package
    % https://github.com/Eomys/MoSQITo/
    % Adapted here as a derivative work under Apache License 2.0
    % https://www.apache.org/licenses/LICENSE-2.0
    
    % TODO: replace the MoSQITo-based loop approach with a parallelised
    % approach!

    modFundRate = zeros(n_blocks, length(halfBark));
    modMaxWeight = zeros(size(modHiWeight));
    for zBand = nBands:-1:1
        waitbar(i_step/n_steps, w, "Calculating modulation rates in 53 bands, ",...
                num2str(zBand), " to go...");...
        i_step = i_step + 1;

        for llBlock = n_blocks:-1:1
            if max(modRate(:, llBlock, zBand)) > 0
                modRateForLoop = modRate(modRate(:, llBlock, zBand) > 0,...
                                         llBlock, zBand);

                NPeaks = length(modRateForLoop);
                I_i0 = {};
                E_i0 = double.empty(NPeaks, 0);
                
                for iiPeak = NPeaks:-1:1
                    modRateRatio = round(modRateForLoop/modRateForLoop(iiPeak));
                    [uniqVals, iiADupes, iiCDupes] = unique(modRateRatio);
                    countDupes = accumarray(iiCDupes, 1);

                    candidateInds = iiADupes(countDupes==1);

                    if max(countDupes) > 1
                        b = uniqVals(countDupes > 1);
                        for jj = length(uniqVals(countDupes > 1)):-1:1
                            c = b(jj);
                            
                            ic = find(modRateRatio == c);
                            crit = abs(modRateForLoop(ic)./(modRateRatio(ic)*modRateForLoop(iiPeak)) - 1);

                            [~, argMin] = min(crit);
                            candidateInds(end + 1) = ic(argMin);
                            candidateInds = fliplr(candidateInds);

                        end
                    end

                    hComplex = abs(modRateForLoop(candidateInds)./(modRateRatio(candidateInds)*modRateForLoop(iiPeak)) - 1);
                    I_i0{iiPeak}= candidateInds(hComplex < 0.04);

                    E_i0(iiPeak) = sum(modHiWeight(I_i0{iiPeak}, llBlock, zBand));
        
                end
            
                [~, i_max] = max(E_i0);
                I_max = I_i0{i_max};
                modFundRate(llBlock, zBand) = modRateForLoop(i_max);
                [~, i_peak] = max(modHiWeight(I_max, llBlock, zBand));
    
                gravityWeight = 1 + 0.1*abs(sum(modRateForLoop(I_max).*modHiWeight(I_max, llBlock, zBand))/sum(modHiWeight(I_max, llBlock, zBand))...
                                            - modRateForLoop(i_peak)).^0.749;
                modMaxWeight(I_max, llBlock, zBand) = modHiWeight(I_max, llBlock, zBand).*gravityWeight;
            
            end
        end
    end

    % end of section attributable to MoSQITo

    % Equation 95
    roughLoWeight = hmSRoughWeight(modFundRate, modfreqMaxWeight, roughLoWeightParams);
    modMaxLoWeight = sum(permute(repmat(roughLoWeight, 1, 1, 10), [3, 1, 2]).*modMaxWeight, 1);
    modMaxWeightSum = sum(modMaxWeight, 1);
    mask = modFundRate >= modfreqMaxWeight;
    modMaxLoWeight(mask) = modMaxWeightSum(mask);
    modMaxLoWeight = squeeze(modMaxLoWeight);
    modAmpMax = modMaxLoWeight;
    modAmpMax(modAmpMax < 0.074376) = 0;

    % Time-dependent specific roughness
    % ---------------------------------
    % Section 7.1.7 ECMA-418-2:2022

    % interpolation to 50 Hz sampling rate
    % Section 7.1.7 Equation 103
    l_50 = ceil(n_samples/sampleRate48k*sampleRate50);
    x = (iBlocks - 1)/sampleRatein;
    xq = linspace(0, signalT, l_50);
    for zBand = nBands:-1:1
        specRoughEst(:, zBand) = pchip(x, modAmpMax(:, zBand), xq);
    end  % end of for loop for interpolation
    specRoughEst(specRoughEst < 0) = 0;

    % Section 7.1.7 Equation 107
    specRoughEstRMS = rms(specRoughEst, 2);

    % Section 7.1.7 Equation 108
    specRoughEstAvg = mean(specRoughEst, 2);

    % Section 7.1.7 Equation 106
    Bl50 = zeros(size(specRoughEstAvg));
    mask = specRoughEstAvg ~= 0;
    Bl50(mask) = specRoughEstRMS(mask)./specRoughEstAvg(mask);

    % Section 7.1.7 Equation 105
    El50 = (0.95555 - 0.58449)*(tanh(1.6407*(Bl50 - 2.5804)) + 1)*0.5 + 0.58449;

    % Section 7.1.7 Equation 104
    specRoughEstTform = 0.0180685*c_Rx*(specRoughEst.^El50); %c_R*c_Rx*specRoughEst.^El50;

    % Section 7.1.7 Equation 110
    timeConstants = 0.5*ones(size(specRoughEstTform));
    mask = diff(specRoughEstTform, 1, 1) >= 0;
    timeConstants(mask) = 0.0625;

    % Section 7.1.7 Equation 109
    specRoughness(:, :, chan) = specRoughEstTform;
    specRoughness(2:end, :, chan) = specRoughEstTform(2:end, :).*(1 - exp(-1./(sampleRate50*timeConstants(2:end, :))))...
                                    + specRoughEstTform(1:end - 1, :).*exp(-1./(sampleRate50*timeConstants(2:end, :)));

    close(w)  % close waitbar

end  % end of for loop over channels

% Binaural roughness
% Section 7.1.11 ECMA-418-2:2022
if inchans == 2 && binaural
    specRoughness(:, :, 3) = sqrt(sum(specRoughness.^2, 3)/2);  % Equation 112
    outchans = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Binaural"];
else
    outchans = inchans;  % assign number of output channels
end

% Section 7.1.8 ECMA-418-2:2022
% Time-averaged specific roughness
specRoughnessAvg = mean(specRoughness(17:end, :, :), 1);

% Section 7.1.9 ECMA-418-2:2022
% Time-dependent roughness Equation 111
% Discard singleton dimensions
if outchans == 1
    roughnessTDep = sum(specRoughness.*dz, 2);
else
    roughnessTDep = squeeze(sum(specRoughness.*dz, 2));
    specRoughnessAvg = squeeze(specRoughnessAvg);
end

% Section 7.1.10 ECMA-418-2:2022
% Overall roughness
roughness90Pc = prctile(roughnessTDep(17:end, :, :), 90, 1);

% time (s) corresponding with results output
timeOut = (0:(size(specRoughness, 1) - 1))/sampleRate50;

%% Output assignment

% Assign outputs to structure
if outchans == 3
    roughnessHMS.specRoughness = specRoughness(:, :, 1:2);
    roughnessHMS.specRoughnessAvg = specRoughnessAvg(:, 1:2);
    roughnessHMS.roughnessTDep = roughnessTDep(:, 1:2);
    roughnessHMS.roughness90Pc = roughness90Pc(:, 1:2);
    roughnessHMS.specRoughnessBin = specRoughness(:, :, 3);
    roughnessHMS.specRoughnessAvgBin = specRoughnessAvg(:, 3);
    roughnessHMS.roughnessTDepBin = roughnessTDep(:, 3);
    roughnessHMS.roughness90PcBin = roughness90Pc(:, 3);
    roughnessHMS.bandCentreFreqs = bandCentreFreqs;
    roughnessHMS.timeOut = timeOut;
else
    roughnessHMS.specRoughness = specRoughness;
    roughnessHMS.specRoughnessAvg = specRoughnessAvg;
    roughnessHMS.roughnessTDep = roughnessTDep;
    roughnessHMS.roughness90Pc = roughness90Pc;
    roughnessHMS.bandCentreFreqs = bandCentreFreqs;
    roughnessHMS.timeOut = timeOut;
end

%% Output plotting

if outplot
    % Plot figures
    % ------------
    for chan = outchans:-1:1
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
        set(get(h,'label'),'string', {'Specific roughness,'; 'asper_{HMS}/Bark_{HMS}'});        
        chan_lab = chans(chan);

        % Create A-weighting filter
        weightFilt = weightingFilter('A-weighting', sampleRatein);
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
             cmap_inferno(34, :), 'LineWidth', 0.75, 'DisplayName', "90th" + string(newline) + "percentile");
        hold on
        plot(ax2, timeOut, roughnessTDep(:, chan), 'color', cmap_inferno(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        if max(roughnessTDep(:, chan)) > 0
            ax2.YLim = [0, 1.05*ceil(max(roughnessTDep(:, chan))*10)/10];
        end
        ax2.XLabel.String = 'Time, s';
        ax2.YLabel.String = 'Roughness, asper_{HMS}';
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.FontName = 'Arial';
        ax2.FontSize = 12;
        lgd = legend('Location', 'eastoutside', 'FontSize', 8);
        lgd.Title.String = "Overall";
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting if outplot true

% end of function
