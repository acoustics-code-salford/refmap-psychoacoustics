function z = hertz2Bark(hz, method)
% z = hertz2Bark(hz, method)
%
%
% Inputs
% ------
% hz : vector
%   The input vector of frequencies.
%
% Returns
% -------
% z : vector
%   The critical band rate vector corresponding with the input frequencies.
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
        hz (1, :) double {mustBeInRange(hz, 0, 20000)}
        method (1, :) string {mustBeMember(method,...
                                           {'volk',...
                                            'traunmueller',...
                                            'schroeder', ...
                                            'zwicker', ...
                                            'sottek'})} = 'volk'
    end


%% Calculation

switch method
    case 'volk'
        z = 32.12*(1 - (1 + (hz/873.47).^1.18).^-0.4);

    case 'traunmueller'
        z = 26.81*hz./(1960 + hz) - 0.53;

        z(z < 2) = z(z < 2) + 0.15.*(2 - z(z < 2));
        z(z > 20.1) = z(z > 20.1) + 0.22.*(z(z > 20.1)- 20.1);

    case 'schroeder'
        z = 7*log(hz/650 + ((hz/650).^2 + 1).^0.5);
    
    case 'zwicker'
        z = 13*atan(0.00076*hz) + 3.5*atan((hz/7500).^2);

    case 'sottek'
        deltaFreq0 = 81.9289;
        c = 0.1618;

        z = 1/c.*asinh(c/deltaFreq0*hz);

end

end

