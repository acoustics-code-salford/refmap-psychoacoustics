function bark = shmHz2Bark(hz)
% bark = shmHz2Bark(hz)
%
% Returns the critical band rate(s) in Bark corresponding with the
% frequency value(s) in Hz, according to the Sottek Hearing Model as
% defined in ECMA-418-2:2025.
%
% Inputs
% ------
% hz : float or array of floats
%   Frequency value(s), corresponding with input bark.
%
% Returns
% -------
% 
% bark: float or array of floats
%   Critical band rate(s).
%
% Assumptions
% -----------
% The output bark is returned based on the inverse of the function defined in ECMA-418-2:2025.
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
% Date created: 26/08/2026
% Date last modified: 26/08/2026
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
        hz (:, :) double {mustBePositive}
    end

%% Define constants
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2025 [deltaf(f=0)]
c = 0.1618;  % frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2025

bark = asinh(hz./((deltaFreq0/c)))./c;  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025, inverted

end
