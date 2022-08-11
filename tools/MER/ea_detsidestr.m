function [sidestr, side_ix, track] = ea_detsidestr(str)

% Small function to parse side inputs
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel

track = '';

str_parts = ea_strsplit(str, '_');
if length(str_parts) == 2
    str = str_parts{2};
    if startsWith(str_parts{1}, 'keycheck', 'IgnoreCase', true)
        track = regexprep(str_parts{1}, '^keycheck', '');
    elseif startsWith(str_parts{1}, 'toggle', 'IgnoreCase', true)
        track = regexprep(str_parts{1}, '^toggle', '');
    elseif startsWith(str_parts{1}, {'key', 'pos'}, 'IgnoreCase', true)
        track = regexprep(str_parts{1}, '^(key|pos))', '');
    end
end

side_ix = [];
sidestr={};
if strcmpi(str,'right')
    side_ix = 1;
    sidestr = {'right',''};
elseif strcmpi(str,'left')
    side_ix = 2;
    sidestr = {'','left'};
elseif strcmpi(str,'both')
    side_ix = 1:2;
    sidestr = {'right','left'};
end
