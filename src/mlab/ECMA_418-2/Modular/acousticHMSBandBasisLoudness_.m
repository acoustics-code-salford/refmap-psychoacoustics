function [signalRectSeg, bandBasisLoudness, bandRMS]...
    = acousticHMSBandBasisLoudness_(signalSegmented, bandCentreFreq)
% [signalRectSeg, bandBasisLoudness, bandRMS]...
% = acousticHMSBandBasisLoudness_(signalSegmented, bandCentreFreq)
%
% Returns rectified input and basis loudness in specified half-Bark
% critical band according to ECMA-418-2:2022 (the Hearing Model of Sottek)
% for an input band-limited signal, segmented into processing blocks
%
% Inputs
% ------
% signalSegmented : 2D matrix
%                   input band-limited segmented signal
%
% bandCentreFreq : double
%                  half-Bark critical band centre frequency
% 
% Returns
% -------
% signalRectSeg : 2D matrix
%                 rectified band-limited segmented signal
%
% bandBasisLoudness : 2D matrix
%                     basis loudness in each block, orientated as per the
%                     input
% 
% bandRMS : column vector
%           RMS for each block (orientated as per the input)
%
% Assumptions
% -----------
% The input signal is a band-limited, segmented signal obtained using
% acousticHMSAuditoryFiltBank.m and signalSegment.m
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (matt@anv.co.uk)
% Institution: University of Salford / ANV Measurement Systems
%
% Date created: 27/09/2023
% Date last modified: 27/09/2023
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
        signalSegmented (:, :) double {mustBeReal}
        bandCentreFreq (1, 1) double {mustBePositive}
    end

%% Define constants

deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2022
c = 0.1618;  % Half-Bark band centre-frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2022

halfBark = 0.5:0.5:26.5;  % half-critical band rate scale
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2022

cal_N = 0.0211668;  % Calibration factor from Section 5.1.8 Equation 23 ECMA-418-2:2022
cal_Nx = 1.0023565;  % Adjustment to calibration factor cal_N Footnote 9 ECMA-418-2:2022
a = 1.5;  % Constant (alpha) from Section 5.1.8 Equation 23 ECMA-418-2:2022

% Values from Section 5.1.8 Table 2 ECMA-418-2:2022
p_threshold = 2e-5*10.^((15:10:85)/20).';
v = [1, 0.6602, 0.0864, 0.6384, 0.0328, 0.4068, 0.2082, 0.3994, 0.6434];

% Loudness threshold in quiet Section 5.1.9 Table 3 ECMA-418-2:2022
LTQz = [0.3310, 0.1625, 0.1051, 0.0757, 0.0576, 0.0453, 0.0365, 0.0298,...
        0.0247, 0.0207, 0.0176, 0.0151, 0.0131, 0.0115, 0.0103, 0.0093,...
        0.0086, 0.0081, 0.0077, 0.0074, 0.0073, 0.0072, 0.0071, 0.0072,...
        0.0073, 0.0074, 0.0076, 0.0079, 0.0082, 0.0086, 0.0092, 0.0100,...
        0.0109, 0.0122, 0.0138, 0.0157, 0.0172, 0.0180, 0.0180, 0.0177,...
        0.0176, 0.0177, 0.0182, 0.0190, 0.0202, 0.0217, 0.0237, 0.0263,...
        0.0296, 0.0339, 0.0398, 0.0485, 0.0622];

%% Input check

if ~ismember(bandCentreFreq, bandCentreFreqs)
    error("Input half-Bark critical rate scale band centre frequency does not match ECMA-418-2:2022 values")
end

%% Signal processing

% Half Wave Rectification
% -----------------------
% Section 5.1.6 Equation 21 ECMA-418-2:2020
signalRectSeg = signalSegmented;
signalRectSeg(signalSegmented <= 0) = 0;

% Calculation of RMS
% ------------------
% Section 5.1.7 Equation 22 ECMA-418-2:2022
bandRMS = sqrt((2/size(signalRectSeg, 1))*sum(signalRectSeg.^2, 1));

% Transformation into Loudness
% ----------------------------
% Section 5.1.8 Equations 23 & 24 ECMA-418-2:2022
bandLoudness = cal_N*cal_Nx*(bandRMS/20e-6).*prod((1 + (bandRMS./p_threshold).^a).^((diff(v)/a)'));
% Section 5.1.9 Equation 25 ECMA-418-2:2022
zBand = find(bandCentreFreq == bandCentreFreqs);
bandBasisLoudness = bandLoudness - LTQz(zBand);  % half-Bark critical band basis loudness
bandBasisLoudness((bandLoudness - LTQz(zBand)) < 0) = 0;

% end of function