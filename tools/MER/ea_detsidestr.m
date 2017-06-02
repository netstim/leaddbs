function [sidestr,side,track] = ea_detsidestr(str)

% Small function to parse side inputs
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel

sidestr={};
side=[];
if strcmpi(str,'right')
    side = 1;
    sidestr = {'right',''};
elseif strcmpi(str,'left')
    side = 2;
    sidestr = {'','left'};
elseif strcmpi(str,'both')
    side=1:2;
    sidestr = {'right','left'};    
end

% Keymer
if isempty(sidestr) || isempty(side) && strcmpi(str(1:3),'key')
    track=str(4:strfind(str,'_')-1);
    tmp = regexp(str,{'right','left'},'match');
    if ~isempty(tmp{1})
        side=1;
        sidestr = {'right',''};
    elseif ~isempty(tmp{2})
        side = 2;
        sidestr = {'','left'};
    end
end