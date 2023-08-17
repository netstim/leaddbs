function atlassurfs = ea_meshatlaslabels(varargin)
%
% Small function to show atlases in mesh form given string


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
        atlassurfs{i}.edgecolor=[1 1 1];
        atlassurfs{i}.color='none';
    end
end

atlassurfs = atlassurfs(idx);
