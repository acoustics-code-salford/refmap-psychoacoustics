function [startIdx, endIdx] = findSequenceIndices(A, mode)
% FINDSEQUENCEINDICES identifies 1s in a 3D Boolean matrix.
%
% [startIdx, endIdx] = findSequenceIndices(A, mode)
%
% Inputs:
%   A    - 3D logical matrix of size [rows, cols, pages]
%   mode - 'longest'  : Find longest sequence of 1s along dim 1
%          'firstlast': Find first and last 1 in each column-plane
%
% Outputs:
%   startIdx - Matrix of start indices (cols x pages)
%   endIdx   - Matrix of end indices   (cols x pages)

[rows, cols, pages] = size(A);
startIdx = NaN(cols, pages);
endIdx   = NaN(cols, pages);

switch lower(mode)
    case 'longest'
        % Pad A along row dimension
        Apad = padarray(A, [1 0 0], 0, 'both');
        starts = diff(Apad, 1, 1) == 1;
        ends   = diff(Apad, 1, 1) == -1;

        % Get start/end indices
        [start_i, start_j, start_k] = ind2sub(size(starts), find(starts));
        [end_i,   end_j,   end_k  ] = ind2sub(size(ends),   find(ends));
        lengths = end_i - start_i;

        % Group by column-plane
        groupID = sub2ind([cols, pages], start_j, start_k);
        T = table(groupID, start_i, end_i, lengths);

        % Sort by group and descending length
        T = sortrows(T, {'groupID', 'lengths'}, {'ascend', 'descend'});
        [uniqueGID, ia] = unique(T.groupID, 'stable');
        Tmax = T(ia, :);

        % Map back to col/page
        [col, page] = ind2sub([cols, pages], Tmax.groupID);
        linIdx = sub2ind([cols, pages], col, page);
        startIdx(linIdx) = Tmax.start_i;
        endIdx(linIdx)   = Tmax.end_i;

    case 'firstlast'
        % Reshape to 2D: rows x (cols*pages)
        A2 = reshape(A, rows, []);
        firstIdxFlat = NaN(1, size(A2,2));
        lastIdxFlat  = NaN(1, size(A2,2));

        % Find first and last 1 per column (vectorized)
        firstMask = A2;
        firstIdxFlat(any(firstMask,1)) = max( firstMask .* (1:rows)', [], 1 );
        lastIdxFlat(any(firstMask,1))  = max( flipud(firstMask) .* (1:rows)', [], 1 );
        lastIdxFlat(any(firstMask,1))  = rows + 1 - lastIdxFlat(any(firstMask,1));

        % Reshape back to [cols, pages]
        startIdx = reshape(firstIdxFlat, [cols, pages]);
        endIdx   = reshape(lastIdxFlat,  [cols, pages]);

    otherwise
        error('Invalid mode. Use ''longest'' or ''firstlast''.');
end
end