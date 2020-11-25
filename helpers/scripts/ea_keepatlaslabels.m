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
    varargin{1} = 'on';
else
    varargin = lower(varargin);
end

for i = 1:length(atlassurfs)
    atlasTag = regexprep(lower(atlassurfs{i}.Tag), '_(left|right|midline|mixed)$', '');
    if strcmp(varargin{1}, 'on') || any(contains(varargin, atlasTag))
        atlassurfs{i}.Visible = 'on';
        colorbuttons(i).State = 'on';
    else
        atlassurfs{i}.Visible = 'off';
        colorbuttons(i).State = 'off';
    end
end
