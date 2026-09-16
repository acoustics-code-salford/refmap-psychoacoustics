function hz = shmBark2Hz(bark)
% hz = shmBark2Hz(bark)
%
% Returns the frequency value(s) in Hz corresponding with the critical
% band rate(s) in Bark, according to the Sottek Hearing Model as defined in
% ECMA-418-2:2025.
%
% Inputs
% ------
% bark : float or array of floats
%   Critical band rate(s).
%
% Returns
% -------
% 
% hz: float or array of floats
%   Frequency value(s), corresponding with input bark.
%
% Assumptions
% -----------
% The output hz is returned based on the function defined in ECMA-418-2:2025.
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
        bark (:, :) double {mustBePositive}
    end

%% Define constants
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2025 [deltaf(f=0)]
c = 0.1618;  % frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2025

hz = (deltaFreq0/c)*sinh(c*bark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025

end
