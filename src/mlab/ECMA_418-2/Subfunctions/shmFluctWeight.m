function fluctWeight = shmFluctWeight(modRate, bandCentreFreq)
% fluctWeight = shmFluctWeight(modRate, bandCentreFreq)
%
% Returns the fluctuation strength band-pass modulation-rate weighting
% w_lh according to ECMA-418-2:2025 (the Sottek Hearing Model), Section
% 9.1.6, Equation 148.
%
% Inputs
% ------
% modRate : array
%   modulation rate(s) [Hz] at which to evaluate the weighting (values
%   of 0 return a weighting of zero)
%
% bandCentreFreq : double
%   critical band centre frequency F(z) [Hz], used in the carrier-
%   frequency correction applied to the high-modulation-rate branch
%
% Returns
% -------
% fluctWeight : array
%   the weighting values w_lh(modRate), the same size as modRate
%
% Assumptions
% -----------
% bandCentreFreq is a scalar value corresponding with a single critical
% band, consistent with a single call of this function per band.
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 16/09/2026
% Date last modified: 16/09/2026
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
%% Arguments validation
    arguments (Input)
        modRate double {mustBeReal, mustBeNonnegative}
        bandCentreFreq (1, 1) double {mustBePositive}
    end

%% Define constants

% Section 9.1.6 Equation 148
fMax = 4.8659;  % [f_max], Hz: modulation rate of maximum (unity) weighting
q1l = 0.33048;
q2l = 0.85902;
q1h = 0.21792;
q2h = 4.6728;

% carrier-frequency correction factor for the high-modulation-rate branch
freqCorrection = (1 + 0.092623*abs(log2(bandCentreFreq/1000))^1.24)^-1;

%% Signal processing

fluctWeight = zeros(size(modRate));

maskLo = modRate > 0 & modRate <= fMax;
maskHi = modRate > fMax;

fluctWeight(maskLo) = (1./(1 + ((modRate(maskLo)./fMax - fMax./modRate(maskLo)).*q1l).^2)).^q2l;
fluctWeight(maskHi) = freqCorrection.*(1./(1 + ((modRate(maskHi)./fMax - fMax./modRate(maskHi)).*q1h).^2)).^q2h;

% end of function