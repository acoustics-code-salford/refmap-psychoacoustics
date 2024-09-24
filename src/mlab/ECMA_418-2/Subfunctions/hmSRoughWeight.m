function roughWeight = hmSRoughWeight(modRate, modfreqMaxWeight, roughWeightParams)
% roughWeight = hmSRoughWeight(modRate, modfreqMaxWeight, roughWeightParams)
%
%
% Inputs
% ------
%
% Returns
% -------
%
%
% Assumptions
% -----------
% Inputs are in compatible parallelised forms
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
% Date created: 10/07/2024
% Date last modified: 24/09/2024
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
% Equation 85
roughWeight = 1./...
                (1 +...
                 ((modRate./modfreqMaxWeight...
                   - modfreqMaxWeight./modRate)...
                  .*roughWeightParams(1, :, :)).^2).^roughWeightParams(2, :, :);