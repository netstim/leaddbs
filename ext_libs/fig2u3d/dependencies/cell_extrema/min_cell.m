function answer = min_cell(a)
%MIN_CELL   minimum of 2D cell or 2D numeric matrix.
%   MIN_CELL(A) calculates the minimum of the contents of A,
%   provided A is a 2 dimesional numeric matrix or
%   a 2 dimensional cell matrix containing a mix of 
%   2D cell matrices and 2D numeric matrices, which in
%   turn may recursively contain others.
%
%   They should not contain any characters or strings,
%   that is, recursive parsing of the elements should
%   end with 2D numeric matrices, empty numeric matrices
%   or empty cell matrices.
%
%   If A is the empty array [] MIN_CELL returns an empty array.
%
%   See also EXTREMA, MAX_CELL, MIN, ISCELL, ISNUMERIC, ISEMPTY.
%
% File:      min_cell.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2009.07.15 - 2010.02.21
% Language:  MATLAB R2011b
% Purpose:   Minimum of 2D Cell/Numeric Matrix
% Copyright: Ioannis Filippidis, 2009-

if iscell(a)
    if isempty(a)
        answer = [];
    else
        a_mins = zeros(1,numel(a));
        k = 1;
        for i=1:size(a,1)
            for j=1:size(a,2)
                b = a{i, j};
                temp = min_cell(b);
                if (isempty(temp) == 0)
                    a_mins(k) = temp;
                    k = k + 1;
                end
            end
        end
        answer = min(a_mins(1:k-1));
    end
elseif isnumeric(a)
    if isempty(a)
        answer = [];
    else
        n = ndims(a);
        for i=1:n
            a = min(a, [], n+1-i);
        end
        answer = a;
    end
else
    disp('Cannot find min. Neither cell nor numeric matrix!')
end
