function signalFadePad = acousticHMSPreProc_(signal, blockSize, hopSize)
% signalFadePad = acousticHMSPreProc_(signal)
%
% Returns signal with fade-in and zero-padding pre-processing according to
% ECMA-418-2:2022 (the Hearing Model of Sottek) for an input signal.
%
% Inputs
% ------
% signal : vector or 2D matrix
%          the input signal/s
%
% blockSize : integer
%             the maximum signal segmentation block size
%
% hopSize : integer
%           the maximum signal segmentation hop size
%           = (1 - overlap)*blockSize
% 
% Returns
% -------
% signalFadePad : vector or 2D matrix
%                 the output faded, padded signal
%
% Assumptions
% -----------
% The input signal is oriented with time on axis 1 (and channel # on axis
% 2), ie, the fade and padding operation is applied along axis 1.
% The input signal must be sampled at 48 kHz.
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 26/09/2023
% Date last modified: 19/10/2023
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
        signal (:, :) double {mustBeReal}
        blockSize (1, 1) {mustBeInteger}
        hopSize (1, 1) {mustBeInteger}
    end


%% Signal processing

% Input pre-processing
% --------------------
%
% Fade in weighting function Section 5.1.2 ECMA-418-2:2022
fadeWeight = repmat(transpose(0.5 - 0.5*cos(pi*(0:239)/240)), 1, size(signal, 2));
% Apply fade in
signalFade = [fadeWeight.*signal(1:240, :);
         signal(241:end, :)];

% Zero-padding Section 5.1.2 ECMA-418-2:2022
n_zeross = max(blockSize);  % start zero-padding
n_samples = size(signal, 1);
n_new = max(hopSize)*(ceil((n_samples + max(hopSize) + n_zeross)/max(hopSize)) - 1);
n_zerose = n_new - n_samples;  % end zero-padding
% Apply zero-padding
signalFadePad = [zeros(n_zeross, size(signalFade, 2));
      signalFade;
      zeros(n_zerose, size(signalFade, 2))];

% end of function