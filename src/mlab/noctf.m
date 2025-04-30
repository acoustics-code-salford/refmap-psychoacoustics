function [fm, f1, f2, fn] = noctf(frange, b)
% Return consecutive range of exact and nominal mid-frequencies, with lower
% and upper band-edge frequencies for 1/b fractional-octave-band filters,
% according to BS EN 61260:2014
% 
% Inputs
% ------
% frange: vector (default = [20 20e3])
%         two numerical values defining the frequency range of interest
% b : integer
%     number defining the octave fraction 1/b
% 
% Returns
% -------
% fm : vector
%      the band mid frequencies
% f1 : vector
%      the lower band-edge frequencies
% f2 : vector
%      the upper band-edge frequencies
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
% Date created: 31/03/2025
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
        frange (1, 2) double {mustBeReal} = [20, 20e3]
        b (1, 1) {mustBeInteger, mustBeGreaterThanOrEqual(b, 1)} = 1
    end
%% Define constants

fref = 1e3;  % reference frequency
octRatio10 = 10^(3/10);  % octave ratio (base-10)

%% Calculate exact frequencies

ind = -15*b:1:15*b;  % large range of indices to calculate frequencies

if mod(b, 1) == 0
    f = octRatio10.^(ind./b)*fref;
else
    f = octRatio10.^((2*ind + 1)./(2*b))*fref;
end

% find range of relevance
[~, ilow] = min(abs(f - min(frange)));  % index of highest f in frange
[~, ihigh] = min(abs(f - max(frange)));  % index of lowest f in frange

% band mid frequencies
fm = f(ilow:ihigh);

% band-edge frequencies
f1 = fm.*octRatio10^(-1/(2*b));
f2 = fm.*octRatio10^(1/(2*b));

%% Calculate nominal mid frequencies
fbase = [100, 125, 160, 200, 250, 315, 400, 500, 630, 800];
fbaseWide = fbase.*10.^ind.';
fbaseWide = sort(reshape(fbaseWide, 1, numel(fbaseWide)));

% find range of relevance
[~, jlow] = min(abs(fbaseWide - min(fm)));  % index of highest f in frange
[~, jhigh] = min(abs(fbaseWide - max(fm)));  % index of lowest f in frange

fn = fbaseWide(jlow:jhigh);
