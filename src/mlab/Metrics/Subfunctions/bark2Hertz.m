function hz = bark2Hertz(z, method)
% hz = bark2Hertz(z, method)
%
%
% Inputs
% ------
% z : vector
%   The input vector of critical band rates.
%
% Returns
% -------
% hz : vector
%   The frequency vector corresponding with the input critical band rates.
%
% Assumptions
% -----------
%
% 
%
% References
% ----------
%
%
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
% Date created: 05/03/2026
% Date last modified: 05/03/2026
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
%% Arguments validation;
    arguments (Input)
        z (1, :) double {mustBeInRange(z, 0, 25)}
        method (1, :) string {mustBeMember(method,...
                                           {'volk',...
                                            'traunmueller',...
                                            'schroeder', ...
                                            'sottek'})} = 'volk'
    end


%% Calculation

switch method
    case 'volk'
        hz = 873.47*((32.12./(32.12 - z)).^2.5 - 1).^(1/1.18);

    case 'traunmueller'
        z(z < 2) = (z(z < 2) - 0.3)./z(z < 2);
        z(z > 20.1) = (z(z > 20.1) + 4.422)./1.22;

        hz = 1960*(z + 0.53)./(26.28 - z);

    case 'schroeder'
        hz = 650.*sinh(z/7);

    case 'sottek'
        deltaFreq0 = 81.9289;
        c = 0.1618;

        hz = (deltaFreq0/c)*sinh(c*z); 

end

end

