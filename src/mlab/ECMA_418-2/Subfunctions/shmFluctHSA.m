function out = shmFluctHSA(envMagSqSpectrBandBlock, modSpecCriterionBankBlock, sampleRate)
% out = shmFluctHSA(envMagSqSpectrBandBlock, modSpecCriterionBankBlock, sampleRate)
% envMagSqSpectrBandBlock : <Description placeholder>
%   <Description placeholder>
%
% modSpecCriterionBankBlock: <Description placeholder>
%   <Description placeholder>
%
% sampleRate : <Description placeholder>
%   <Description placeholder>
%
% Returns
% -------
% out : <Description placeholder>
%   <Description placeholder>
%
% Assumptions
% -----------
% <Description placeholder>
%
% Requirements
% ------------
% <Description placeholder>
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 31/10/2025
% Date last modified: 02/11/2025
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
% <Description placeholder>

mlabIndex = 1;
% identify peaks in each block (for each band)
% NOTE: peak search limited to k = 0,...,48
[PhiPks, pkLocs, ~, ~] = findpeaks(envMagSqSpectraBandBlock(0 + mlabIndex:48 + mlabIndex),...
                                  'MinPeakHeight', modSpecCriterionBandBlock);

if ~isempty(PhiPks)

    kLocs = pkLocs - mlabIndex;  % 0-indexed k values
    % modRate(length(PhiPksMask),...
    %         lBlock,...
    %         zBand) = 

    deltaF1500*

end  % end of if branch for thresholded modulation spectra


