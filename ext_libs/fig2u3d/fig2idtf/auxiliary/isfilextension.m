function [yes] = isfilextension(fname, extension)
% See also CHECK_FILE_EXTENSION, CLEAR_FILE_EXTENSION.
%
% File:      isfilextension.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.22
% Language:  MATLAB R2012a
% Purpose:   check if file has the desired extension
% Copyright: Ioannis Filippidis, 2012-

% dot in extension ?
if ~strcmp(extension(1), '.')
    extension = ['.', extension];
end

[~, ~, ext] = fileparts(fname);

% extension = unwanted ?
if strcmp(ext, extension)
    yes = 1;
else
    yes = 0;
end
