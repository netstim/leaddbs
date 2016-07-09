function answer = max_cell(a)
%MAX_CELL   maximum of 2D cell or 2D numeric matrix.
%   MAX_CELL(A) calculates the maximum of the contents of A,
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
%   If A is the empty array [] MAX_CELL returns an empty array.
%
%   See also EXTREMA, MIN_CELL, MAX, ISCELL, ISNUMERIC, ISEMPTY.
%
% File:      max_cell.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2009.07.15 - 2010.02.21
% Language:  MATLAB R2011b
% Purpose:   Maximum of 2D Cell/Numeric Matrix
% Copyright: Ioannis Filippidis, 2009-

if iscell(a)
    if isempty(a)
        answer = [];
    else
        a_maxs = zeros(1,numel(a));
        k = 1;
        for i=1:size(a,1)
            for j=1:size(a,2)
                b = a{i, j};
                temp = max_cell(b);
                if (isempty(temp) == 0)
                    a_maxs(k) = temp;
                    k = k + 1;
                end
            end
        end
        answer = max(a_maxs(1:k-1));
    end
elseif isnumeric(a)
    if isempty(a)
        answer = [];
    else
        n = ndims(a);
        for i=1:n
            a = max(a, [], n+1-i);
        end
        answer = a;
    end
else
    disp('Cannot find max. Neither cell nor numeric matrix!')
end
