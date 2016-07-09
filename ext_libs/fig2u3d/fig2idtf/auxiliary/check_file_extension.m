function [fname] = check_file_extension(fname, extension)
% See also CLEAR_FILE_EXTENSION, ISFILEXTENSION.
%
% File:      check_file_extension.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.22
% Language:  MATLAB R2012a
% Purpose:   append desired file extension if absent
% Copyright: Ioannis Filippidis, 2012-

% extension ok ?
if isfilextension(fname, extension)
    return
end

% dot in extension ?
if ~strcmp(extension(1), '.')
    extension = ['.', extension];
end

fname = [fname, extension];
