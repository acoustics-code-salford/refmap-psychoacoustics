function soundLevel = loudLfISO226_2023(alphaf, magTransfer, thresholdf, phonLevel, axisN)
% soundLevel = loudLfISO226_2023(loudExponent, magTransfer, threshold)
% 
% Return sound pressure level(s) of a pure tone according to ISO 226:2023
% (for 1/3 octave bands)
%
% Inputs
% ------
% alphaf : double or vector
%          exponent(s) of loudness perception over frequencies.
%
% magTransfer : double or vector
%               normalised magnitude of linear transfer function (L_U).
%
% thresholdf : double or vector
%              hearing threshold (Tf).
%
% phonLevel : double, 
%             input loudness level
%
% axisN : integer (default: 1)
%         series axis (that is not frequency bands)
%
% Returns
% -------
% soundLevel : double or vector
%              sound pressure level(s)
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
% Date created: 22/04/2025
% Date last modified: 22/04/2025
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
        alphaf (1, :) double {mustBeReal}
        magTransfer (1, :) double {mustBeReal}
        thresholdf (1, :) double {mustBeReal}
        phonLevel (:, :) double {mustBeReal}
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
    end

%% Input check
if axisN == 2
    alphaf = alphaf.';
    magTransfer = magTransfer.';
    thresholdf = thresholdf.';
    phonLevel = phonLevel.';
    axisTpose = true;
else
    axisTpose = false;
end  % end of input orientation if branch

% if more than one loudness level, check length matches other inputs
if size(phonLevel, 2) > 1
    if ~isequal(size(phonLevel, 2), size(alphaf, 2), size(magTransfer, 2), size(thresholdf, 2))
        error("Input variables must have the same number of frequency bands")
    end
% otherwise, check remaining inputs have same length
else
    if ~isequal(size(alphaf, 2), size(magTransfer, 2), size(thresholdf, 2))
        error("Input variables must have the same number of frequency bands")
    end
end

% if a series of input loudness levels
if size(phonLevel, 1) > 1
    % repeat matrices
    nSpectra = size(phonLevel, 1);
    alphaf = repmat(alphaf, nSpectra, 1);
    magTransfer = repmat(magTransfer, nSpectra, 1);
    thresholdf = repmat(thresholdf, nSpectra, 1);
end

%% Calculation
% ISO 226:2023 Equation 1
soundLevel = 10./alphaf.*log10(4e-10.^(0.3 - alphaf).*(10.^(0.03*phonLevel)...
                                                       - 10^0.072) + 10.^(alphaf.*(thresholdf + magTransfer)/10)) - magTransfer;

% if inputs transposed, reorientate output
if axisTpose
    soundLevel = soundLevel.';
end

end

