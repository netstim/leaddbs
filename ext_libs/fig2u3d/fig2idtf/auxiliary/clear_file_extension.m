function [fname] = clear_file_extension(fname, unwanted_ext)
% See also CHECK_FILE_EXTENSION, ISFILEXTENSION.
%
% File:      clear_file_extension.m
% Author:    Ioannis Filippidis, jfilippidis@gmail.com
% Date:      2012.06.22
% Language:  MATLAB R2012a
% Purpose:   remove the unwanted file name extension, if present
% Copyright: Ioannis Filippidis, 2012-

% extension = unwanted ?
if ~isfilextension(fname, unwanted_ext)
    return
end

%% clear unwanted extension
[path, name] = fileparts(fname);

% preceding path ?
if isempty(path)
    fname = name;
else
    fname = [path, filesep, name];
end
