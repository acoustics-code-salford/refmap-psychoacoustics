function [start_idx, end_idx] = findColumnOneIndices(A, mode)
% FINDCOLUMNONEINDICES Identifies indices of 1s in each column of a 2D/3D logical matrix.
%
%   [start_idx, end_idx] = findColumnOneIndices(A)
%   [start_idx, end_idx] = findColumnOneIndices(A, mode)
%
%   A     : 2D or 3D logical matrix
%   mode  : 'firstlast' (default) - first and last occurrence of 1s
%           'longest'             - longest continuous sequence of 1s
%
%   start_idx, end_idx : (n x p) matrices of row indices for each (col, page)
%                        pair in the input matrix A

    if nargin < 2
        mode = 'firstlast';
    end

    if ~islogical(A)
        error('Input matrix A must be logical.');
    end

    dims = ndims(A);

    if dims == 2
        [m, n] = size(A);
        p = 1;
    elseif dims == 3
        [m, n, p] = size(A);
    else
        error('Input matrix A must be 2D or 3D.');
    end

    total_cols = n * p;

    switch lower(mode)
        case 'firstlast'
            [i, j, k] = ind2sub(size(A), find(A));
            if dims == 2
                k = ones(size(j));
            end
            col_plane_idx = sub2ind([n, p], j, k);
            start_flat = accumarray(col_plane_idx, i, [total_cols, 1], @min, NaN);
            end_flat   = accumarray(col_plane_idx, i, [total_cols, 1], @max, NaN);

        case 'longest'
            % Reshape A to [m, n*p] for easier vectorized processing
            A_reshaped = reshape(A, m, []);
            % Pad with zeros at top and bottom to detect edges
            padded = [zeros(1, total_cols); A_reshaped; zeros(1, total_cols)];
            % Find starts and ends of sequences
            diffed = diff(padded);
            starts = find(diffed == 1);
            ends   = find(diffed == -1) - 1;

            col_idx = mod(starts-1, total_cols) + 1;
            run_len = ends - starts + 1;

            % Find longest run per column
            max_run = accumarray(col_idx, run_len, [total_cols, 1], @max, NaN);
            start_flat = NaN(total_cols,1);
            end_flat   = NaN(total_cols,1);

            % For each column, find the first run matching max length
            valid = ~isnan(max_run);
            for kcol = find(valid)'
                matches = (col_idx == kcol) & (run_len == max_run(kcol));
                first_match = find(matches, 1, 'first');
                start_flat(kcol) = starts(first_match);
                end_flat(kcol)   = ends(first_match);
            end

        otherwise
            error('Invalid mode. Use ''firstlast'' or ''longest''.');
    end

    % Reshape to n x p
    start_idx = reshape(start_flat, [n, p]);
    end_idx   = reshape(end_flat, [n, p]);

    if dims == 2
        start_idx = start_idx(:,1);
        end_idx   = end_idx(:,1);
    end
end