function [sidestr, side_ix, track] = ea_detsidestr(str)

% Small function to parse side inputs
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel

track = '';

str_parts = split(str, '_');
if length(str_parts) == 2
    str = str_parts{2};
    if (length(str_parts{1}) > 8) && strcmpi(str_parts{1}(1:8), 'keycheck')
        track = str_parts{1}(9:end);
    elseif (length(str_parts{1}) > 6) && strcmpi(str_parts{1}(1:6), 'toggle')
        track = str_parts{1}(7:end);
    elseif (length(str_parts{1}) > 3) && any(strcmpi(str_parts{1}(1:3), {'key', 'pos'}))
        track = str_parts{1}(4:end);
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