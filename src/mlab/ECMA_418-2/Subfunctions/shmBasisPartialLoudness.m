function [signalRectSegTarget, basisPartLoudness, blockRMSTarget] = shmBasisPartialLoudness(signalSegmentedTarget, signalSegmentedMasker, bandCentreFreq)
% [signalRectSeg, basisPartLoudness, blockRMS] = shmBasisPartialLoudness(signalSegmentedTarget,
%                                                                        signalSegmentedMasker
%                                                                        bandCentreFreq)
%
% Returns rectified input target signal and basis partial loudness in specified
% critical band based on ECMA-418-2:2025 (the Sottek Hearing Model)
% for an input band-limited target and masker signals, segmented into processing blocks
%
% Inputs
% ------
% signalSegmentedTarget : 2D or 3D matrix
%   Input band-limited segmented target signal(s)
%
% signalSegmentedMasker : 2D or 3D matrix
%   Input band-limited segmented masker signal(s)
%
% bandCentreFreq : double (optional, default = [])
%   critical band centre frequency - if empty, all bands are assumed to be
%   present in the input segmented signal matrix
% 
% Returns
% -------
% signalRectSeg : 2D or 3D matrix
%   rectified band-limited segmented signal
%
% basisPartLoudness : 2D or 3D matrix
%   basis partial loudness in each block
% 
% blockRMS : column vector or 2D matrix
%   RMS for each block
%
% Assumptions
% -----------
% The input signal is a segmented signal (either band-limited, or arranged
% with critical bands over the third dimension) obtained using
% shmAuditoryFiltBank.m and shmSignalSegment.m
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (m.c.torjussen@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 09/12/2023
% Date last modified: 12/12/2025
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
% This code was developed from an original file 'SottekTonality.m' authored
% by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
% The original code has been reused and updated here with permission.
%
% This code includes an expression obtained from the loudness package
% developed by Dominic Ward
% (https://github.com/deeuu/loudness/tree/master), which is re-used under
% the applicable GPL-3.0 licence terms.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        signalSegmentedTarget double {mustBeReal}
        signalSegmentedMasker double {mustBeReal}
        bandCentreFreq double {mustBePositive} = []
    end

% check if input is 2D and includes band centre frequency - otherwise raise
% error
if isempty(bandCentreFreq) && length(size(signalSegmented)) == 2
    error("Band centre frequency must be specified for single band-limited input signal")
end

% check if input band centre frequency is not a vector (arguments
% validation does not allow empty default with specified size)
if ~isempty(bandCentreFreq) && max(size(bandCentreFreq)) ~= 1
    error("Band centre frequency input must be a single value")
end

%% Define constants

deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2025
c = 0.1618;  % Half-Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2025

halfBark = 0.5:0.5:26.5;  % half-critical band rate scale
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025

cal_N = 0.0211668;  % Calibration factor from Section 5.1.8 Equation 23 ECMA-418-2:2025
cal_Nx = 1.00132;  % Calibration multiplier (Footnote 8 ECMA-418-2:2025)

a = 1.5;  % Constant (alpha) from Section 5.1.8 Equation 23 ECMA-418-2:2025

% Values from Section 5.1.8 Table 2 ECMA-418-2:2025
p_threshold = 2e-5*10.^((15:10:85)/20).';
v = [1, 0.6602, 0.0864, 0.6384, 0.0328, 0.4068, 0.2082, 0.3994, 0.6434];

% Loudness threshold in quiet Section 5.1.9 Table 3 ECMA-418-2:2025
LTQz = [0.3310, 0.1625, 0.1051, 0.0757, 0.0576, 0.0453, 0.0365, 0.0298,...
        0.0247, 0.0207, 0.0176, 0.0151, 0.0131, 0.0115, 0.0103, 0.0093,...
        0.0086, 0.0081, 0.0077, 0.0074, 0.0073, 0.0072, 0.0071, 0.0072,...
        0.0073, 0.0074, 0.0076, 0.0079, 0.0082, 0.0086, 0.0092, 0.0100,...
        0.0109, 0.0122, 0.0138, 0.0157, 0.0172, 0.0180, 0.0180, 0.0177,...
        0.0176, 0.0177, 0.0182, 0.0190, 0.0202, 0.0217, 0.0237, 0.0263,...
        0.0296, 0.0339, 0.0398, 0.0485, 0.0622];

% Partial masking parameters
% Approx K from Moore & Glasberg, 1997 (K is the post auditory filter
% signal-to-noise ratio required for threshold excitation in high masking
% noise, here denoted snrTMN)
% Expression adapted from https://github.com/deeuu/loudness/blob/master/src/modules/SpecificPartialLoudnessMGB1997.cpp
% by Dominic Ward
snrTMNdB = ones(size(bandCentreFreqs)).*-3;
mask = bandCentreFreqs < 1000;
logFreq = log10(bandCentreFreqs);
% NOTE: the clamp at bandCentreFreqs <= 50 Hz applied in the original code
% [logFreq(bandCentreFreqs <= 50) = log10(50);] is disregarded here in
% favour of continuing the extrapolation down to the 41 Hz critical band.
snrTMNdB(mask) = 4.4683421395470813*logFreq(mask).^4 ...
                 - 49.680265975480935*logFreq(mask).^3 ...
                 + 212.89091871289099*logFreq(mask).^2 ...
                 - 418.88552055717832*logFreq(mask)...
                 + 317.07466358701885;
snrTMN = 10.^(snrTMNdB/10);
% End of re-used code from loudness package.

% limits for masking function parameters
maskCompressLim = [0.2, 0.25];
maskResponseLim = [1, 2];

% masking function compression parameter
maskCompress = maskCompressLim(1)...
               + (maskCompressLim(2)...
                  - maskCompressLim(1)).*(snrTMN./(1 + snrTMN));

% masking function response parameter
maskResponse = maskResponseLim(1)...
               + (maskResponseLim(2)...
                  - maskResponseLim(1)).*(snrTMN./(1 + snrTMN));

% standardised epsilon
epsilon = 1e-12;

%% Input check

if ~isempty(bandCentreFreq) && ~ismember(bandCentreFreq, bandCentreFreqs)
    error("Input half-Bark critical rate scale band centre frequency does not match ECMA-418-2:2025 values")
end

%% Signal processing

% Half Wave Rectification
% -----------------------
% Section 5.1.6 Equation 21 ECMA-418-2:2020
signalRectSegTarget = signalSegmentedTarget;
signalRectSegTarget(signalSegmentedTarget <= 0) = 0;
signalRectSegMasker = signalSegmentedMasker;
signalRectSegMasker(signalSegmentedMasker <= 0) = 0;

signalSegmentedTotal = signalSegmentedTarget + signalSegmentedMasker;
signalRectSegTotal = signalSegmentedTotal;
signalRectSegTotal(signalSegmentedTotal <= 0) = 0;

% Calculation of RMS
% ------------------
% Section 5.1.7 Equation 22 ECMA-418-2:2025
blockRMSTarget = sqrt((2/size(signalRectSegTarget, 1))*sum(signalRectSegTarget.^2, 1));
blockRMSMasker = sqrt((2/size(signalRectSegMasker, 1))*sum(signalRectSegMasker.^2, 1));
blockRMSTotal = sqrt((2/size(signalRectSegTotal, 1))*sum(signalRectSegTotal.^2, 1));

% Transformation into Loudness
% ----------------------------
% Section 5.1.8 Equations 23 & 24 ECMA-418-2:2025
bandLoudnessTarget = cal_N*cal_Nx*(blockRMSTarget/20e-6).*prod((1 + (blockRMSTarget./p_threshold).^a).^((diff(v)/a)'));
bandLoudnessMasker = cal_N*cal_Nx*(blockRMSMasker/20e-6).*prod((1 + (blockRMSMasker./p_threshold).^a).^((diff(v)/a)'));
bandLoudnessTotal = cal_N*cal_Nx*(blockRMSTotal/20e-6).*prod((1 + (blockRMSTotal./p_threshold).^a).^((diff(v)/a)'));

% remove singleton dimension from outputs
blockRMSTarget = squeeze(blockRMSTarget);
bandLoudnessTarget = squeeze(bandLoudnessTarget);
bandLoudnessMasker = squeeze(bandLoudnessMasker);
bandLoudnessTotal = squeeze(bandLoudnessTotal);

% Section 5.1.9 Equation 25 ECMA-418-2:2025
if ~isempty(bandCentreFreq) && length(size(signalSegmentedTarget)) == 2

    [~, bandIdx] = max(bandCentreFreq == bandCentreFreqs);

    bandLoudnessTarget = max(0, bandLoudnessTarget);
    bandLoudnessMasker = max(0, bandLoudnessMasker);

    % compressed inverse
    invT = bandLoudnessTarget.^(1./maskCompress(bandIdx));
    invM = bandLoudnessMasker.^(1./maskCompress(bandIdx));

    % compressed difference ratio
    compDiffRatio = (invT - invM) ./ (invT + invM + epsilon);
    compDiffRatioPos = max(0, compDiffRatio);
    
    % damping function for compressed difference ratio
    compDiffRatioPosDamped = compDiffRatioPos.^maskResponse(bandIdx);
    
    % partial loudness without LTQ
    bandPartLoudness = (compDiffRatioPosDamped.*invT).^maskCompress(bandIdx);

    % basis partial loudness
    basisPartLoudness = max(0, bandPartLoudness - LTQz(bandIdx));
    basisPartLoudness(bandLoudnessTotal <= LTQz(bandIdx)) = 0;

else

    bandLoudnessTarget = max(0, bandLoudnessTarget);
    bandLoudnessMasker = max(0, bandLoudnessMasker);

    % raise parameter dimensions
    LTQz3D = reshape(repmat(LTQz, size(bandLoudnessMasker, 1)*size(bandLoudnessMasker, 2), 1), size(bandLoudnessMasker));
    maskCompress3D = reshape(repmat(maskCompress, size(bandLoudnessMasker, 1)*size(bandLoudnessMasker, 2), 1), size(bandLoudnessMasker));
    maskResponse3D = reshape(repmat(maskResponse, size(bandLoudnessMasker, 1)*size(bandLoudnessMasker, 2), 1), size(bandLoudnessMasker));

    % compressed inverse
    invT = bandLoudnessTarget.^(1./maskCompress3D);
    invM = bandLoudnessMasker.^(1./maskCompress3D);

    % compressed difference ratio
    compDiffRatio = (invT - invM) ./ (invT + invM + epsilon);
    compDiffRatioPos = max(0, compDiffRatio);
    
    % damping function for compressed difference ratio
    compDiffRatioPosDamped = compDiffRatioPos.^maskResponse3D;
    
    % partial loudness without LTQ
    bandPartLoudness = (compDiffRatioPosDamped.*invT).^maskCompress3D;

    % basis partial loudness
    basisPartLoudness = max(0, bandPartLoudness - LTQz3D);
    basisPartLoudness(bandLoudnessTotal <= LTQz3D) = 0;

end

% end of function
