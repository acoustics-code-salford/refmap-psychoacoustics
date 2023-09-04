function Xwf = fweight_Xf(f, Xf, w_type, axisn)
% [f, Xwf] = fweight_Xf(f, Xf, w_type)
% Returns weighted spectrum from input (unweighted) spectrum.
% (NB: the output X orientation matches the input orientation).
%
% Inputs
% ------
% f : vector of floats
%     the frequencies corresponding with the input spectral values
% Xf : vector or 2D matrix of floats
%      the spectrum or spectra to be weighted (in physical magnitude units,
%      not dB)
% w_type : keyword string
%          the weighting to be applied
% axisn : integer (default 1)
%         the axis along which the spectrum lies (1 or 2)
%
% Outputs
% -------
% Xwf : vector or 2D matrix of floats
%       the weighted spectral values, in the same orientation as Xf
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
% Date created: 22/07/2023
% Date last modified: 13/08/2023
% MATLAB version: 2022b
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
        f (1, :) double {mustBeNonnegative}
        Xf (:, :) double {mustBeNonnegative}
        w_type (1, 1) {mustBeMember(w_type, {'A', 'C'})}
        axisn (1, 1) {mustBeMember(axisn, [1, 2])} = 1
    end

    arguments (Output)
        Xwf (:, :)
    end

%% processing
% check input orientation
if axisn == 2
    transpose(Xf);
end

% apply weighting (as physical unit type, not dB)
switch upper(w_type)
    case 'A'
        weightdB = transpose(20.*log10(12194^2*f.^4./((f.^2 + 20.6^2)...
                             .*sqrt((f.^2 + 107.7^2).*(f.^2 + 737.9^2))...
                             .*(f.^2 + 12194^2))) + 2);
        weight = 10.^(weightdB./20);
        weight(f==0) = 1 - sum(weight(2:end))/sum(ones(size(f(2:end))));  % handling for 0-frequency (DC bias)

    case 'C'
        weightdB = transpose(20.*log10(12194^2*f.^2./((f.^2 + 20.6^2)...
                             .*(f.^2 + 12194^2))) + 0.061904281992313);
        weight = 10.^(weightdB./20);
        weight(f==0) = 1 - sum(weight(2:end))/sum(ones(size(f(2:end))));  % handling for 0-frequency (DC bias)
end

Xwf = Xf.*weight;

% reorientation
if axisn == 2
    transpose(Xwf);
end
