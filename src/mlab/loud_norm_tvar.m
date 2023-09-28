function loud_norm_tvar(ster_type, pcex, iter_lim)
% loud_norm_tvar
% loud_norm_tvar(infilepath, cal_val, cal_type, outfilename)
%
% Saves a set of audio (wav) files calibrated to normalise the total
% loudness of all files, based on the nominal total time-varying loudness
% of the first file in the selection, specified as the percent exceeded of
% the time-varying total loudness (ISO 532-1:2017).
% The normalisation of individual files is arbitrary, since a signal can be 
% reproduced at any desired calibration gain - the normalisation has 
% meaning only when multiple files are processed together and the output 
% files reproduced at the same gain.
% 
% Inputs
% ------
% ster_type : keyword string (optional, default = 'average')
%             the approach to use for equalising the loudness of stereo 
%             signals. Options comprise 'left', 'right', or 'average'
%             (default), indicating which channel is used as the
%             equalisation target.
%
% pcex : integer (optional, 0-100, default = 5)
%        the 'percent exceeded' value for the total loudness calculation
%        (note this is the percent exceeded, and not the statistical
%        percentile value of the cumulative distribution function, which is
%        100 - pcex)
%
% iter_lim : integer (optional, default = 20)
%            the iteration limit on the convergence loop (set to avoid
%            excessive calculation time - increase if needed to converge)
% 
% Returns
% -------
% None
%
% Assumptions
% -----------
% Assumes input and output file formats will be mono or stereo audio in
% .wav format
% The acceptable tolerance for loudness normalisation is assumed to be 0.1
% sones.
%
% Requirements
% ------------
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 24/07/2023
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

%% Arguments validation
    arguments (Input)
        ster_type string {mustBeMember(ster_type, {'left', 'right',...
                          'average'})} = 'average'
        pcex (1, :) double {mustBeNonnegative, mustBeInteger,...
                            mustBeLessThanOrEqual(pcex, 100)} = 95
        iter_lim (1, 1) double {mustBePositive, mustBeInteger} = 20
    end

%% Input section

% open reference file selection dialog
[files1, ~] = uigetfile({'*.wav', "Wave files (*.wav)"},...
                               "Select input audio file to use as normalisation reference");

if isa(files1, 'double')
    disp("Error: No normalisation reference file selected - exiting")
    return
end

% open file selection dialog
[files2, dirpath] = uigetfile({'*.wav', "Wave files (*.wav)"},...
                               "Select input audio files to normalise",...
                               Multiselect='on');

if isa(files2, 'double')
    disp("Error: No files selected for normalisation - exiting")
    return
end

% compile filenames into an array
for ii = length(files2) + 1:-1:1
    if ii > 1
        files(ii) = files2(ii - 1);
    else
        files(ii) = cellstr(files1);
    end
end

% join filenames with directory path
for ii = length(files):-1:1
    filepaths(ii) = join([dirpath, files(ii)], "");
end

%% Processing section
skip_first = false;  % generate truth flag for skipping first file (see below)
for ii = 1:length(files)
    filepath = string(filepaths(ii));
    filename = string(files(ii));

    [x, fs] = audioread(filepath);
    nbits = audioinfo(filepath).BitsPerSample;

    % check input signal channel number and skip if more than two
    if size(x, 2) == 1
        ster_sig = false;
    elseif size(x, 2) == 2
        ster_sig = true;
    else
        fprintf("Warning: Input signals must be mono (1-channel) or stereo (2-channel) - skipping " + filename + "\n");
        if ii == 1
            skip_first = true;
        end
        continue
    end

    % check input signal length and warn if signal is long
    if size(x, 1)/fs > 120
        sig_T = round(size(x, 1)/fs, 0);
        formatSpec = "Warning: Input signal is %i seconds long - expect loudness calculation to take a while!\n";
        fprintf(formatSpec, sig_T);
    end

    % calculate percentile loudness
    [~, ~, percN] = acousticLoudness(x, fs, 1, TimeVarying=true,...
                                     Percentiles=pcex);
    
    % iteration loop to calibrate equal loudness
    % get loudness calibration from first input file and resave output
    if ii == 1 || skip_first
        % for stereo signal, select calibration channel
        if ster_sig
            if strcmp(ster_type, 'left')
                cal_N = percN(1);
            elseif strcmp(ster_type, 'right')
                cal_N = percN(2);
            else
                cal_N = mean(percN);
            end
        else
            cal_N = percN;
        end

        % generate waitbar
        wait_msg = sprintf("Processing %s...", filename);  % waitbar message
        wait_accum = 0;  % waitbar progress initialisation
        h = waitbar(wait_accum, wait_msg);  % waitbar initialisation
        set(findall(h, 'type', 'text'), 'Interpreter', 'none');
        waitbar(wait_accum);  % waitbar update (to fix filename text)

        % write to file
        cal_N_s = num2str(round(cal_N, 1));
        if ster_sig && ~strcmp(ster_type, 'average')
            cal_N_s = strrep(cal_N_s, "         ", "son_");
        end
        outfilename = join([filename, "_", cal_N_s, "son.wav"], "");
        outfilepath = join([dirpath, outfilename], "");
        audiowrite(outfilepath, x, fs, BitsPerSample=nbits);
        
        % switch first file skip truth flag back
        skip_first = false;
    
        % simulated waitbar
        waitbar(0.85);  % update waitbar
        pause(0.5);  % insert pause to enable waitbar to be seen
        waitbar(1);  % update waitbar
        pause(0.25);  % insert pause to enable waitbar to be seen
        close(h);  % close waitbar

    % iterate to adjust signal calibration until loudness matches target to
    % within 0.1 sones.
    else
        % for stereo signal, select normalisation channel
        if ster_sig
            if strcmp(ster_type, 'left')
                percN = percN(1);
            elseif strcmp(ster_type, 'right')
                percN = percN(2);
            else
                percN = mean(percN);
            end
        end

        iters = 1;  % iteration counter
        percN0 = percN;  % store initial loudness value for waitbar
        wait_msg = sprintf("Processing %s...", filename);  % waitbar message
        wait_accum = 0;  % waitbar progress initialisation
        h = waitbar(wait_accum, wait_msg);  % waitbar initialisation
        set(findall(h, 'type', 'text'), 'Interpreter', 'none');
        waitbar(wait_accum);  % waitbar update (to fix filename text)
        while abs(percN - cal_N) > 0.05 && iters <= iter_lim  % iteration condition            
            x = x./(percN/cal_N);  % adjust signal towards target
            [~, ~, percN] = acousticLoudness(x, fs, TimeVarying=true,...
                                       Percentiles=pcex);  % recalculate loudness
            
            % for stereo signal, reselect normalisation channel
            if ster_sig
                if strcmp(ster_type, 'left')
                    percN = percN(1);
                elseif strcmp(ster_type, 'right')
                    percN = percN(2);
                else
                    percN = mean(percN);
                end
            end

            iters = iters + 1;  % increment iteration counter
            wait_accum = abs(percN0 - cal_N) - abs(percN - cal_N);  % update waitbar progress
            waitbar(wait_accum/abs(percN0 - cal_N));  % update waitbar
        end
        if iters > iter_lim  % check if iterations exceeds limit and skip file if so
            formatSpec = "Warning: %s failed to converge within %i iterations - skipping...\n";
            fprintf(formatSpec, filename, iter_lim);
            close(h);
            continue
        end
        waitbar(1);  % update waitbar
        close(h);
        % save adjusted output file
        percN_s = num2str(round(percN, 1));
        percN_s = strrep(percN_s, "         ", "son_");
        outfilename = join([filename, "_", percN_s, "son.wav"], "");
        outfilepath = join([dirpath, outfilename], "");
        audiowrite(outfilepath, x, fs, BitsPerSample=nbits);
        
    end
end
