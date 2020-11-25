function atlassurfs = ea_keepatlaslabels(varargin)
%
% Small function to keep atlas labels given string
%
% Example:
%
%   ea_keepatlaslabels('on');
%   ea_keepatlaslabels('off')
%
%   ea_keepatlaslabels('STN','RN','GPi','GPe')
%   atlassurfs = ea_keepatlaslabels('STN','RN','GP')

% _________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel

H = findall(0,'type','figure');
resultfig = H(contains({H(:).Name},{'Electrode-Scene'}));
resultfig=resultfig(1); % Take the first if there are more.
set(0,'CurrentFigure',resultfig)

atlassurfs = getappdata(resultfig,'atlassurfs');
colorbuttons = getappdata(resultfig,'colorbuttons');

% Show all atlases in case no input parameter
if isempty(varargin)
    atlasToggle = {'on'};
else
    atlasToggle = lower(varargin);
end

idx = [];
for i = 1:length(atlassurfs)
    atlasTag = lower(atlassurfs{i}.Tag);
    atlasName = regexprep(atlasTag, '_(left|right|midline|mixed)$', '');
    if strcmp(atlasToggle{1}, 'on') ... % Turn on all atlases
            || ismember(atlasTag, atlasToggle) ... % Match exactly
            || ismember(atlasName, atlasToggle) % Match name regardless of side
        idx = [idx, i];
        atlassurfs{i}.Visible = 'on';
        colorbuttons(i).State = 'on';
    else
        atlassurfs{i}.Visible = 'off';
        colorbuttons(i).State = 'off';
    end
end

atlassurfs = atlassurfs(idx);
