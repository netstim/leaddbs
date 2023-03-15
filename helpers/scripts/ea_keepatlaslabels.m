function ea_keepatlaslabels(varargin)
%
% Helper function to keep atlas labels or turn them off
%
%
% Example:
%   Default to all objects/atlassurfs in scene
%
%   ea_keepatlaslabels('on');
%   ea_keepatlaslabels('off');
%
% Example 2: include type of object either before or after 'on'/'off'
%    type of object:
%           1. atlas, atlasees, or atlassurfs
%           2. object, objects, roi, rois
%
%   ea_keepatlaslabels('off','rois');
%   ea_keepatlaslabels('off','atlases');
%   ea_keepatlaslabels('rois','on');
%   ea_keepatlaslabels('on','atlas');
% _________________________________________________________________________
% Copyright (C) 2023 Brigham and Women's Hospital
%
% Ari Kappel

% Parse inputs
% Default: Show all if no input parameters provided
try
    if isempty(varargin)
        toggle = 'on';
        type = 'all';
    elseif nargin==1
        toggle = lower(varargin{1});
        type = 'all'; % default to all
    elseif nargin==2
        if any(strcmpi(varargin{1},{'on','off'}))
            toggle = lower(varargin{1});
            type = lower(varargin{2});
        elseif any(strcmpi(varargin{2},{'on','off'}))
            type = lower(varargin{1});
            toggle = lower(varargin{2});
        end
    else
        toggle = 'on';
        type='all';
    end

    % get figure handles
    % Warning: may not work with mutliple scenes open
    H = findall(0,'type','figure');
    resultfig = H(contains({H(:).Name},{'Electrode-Scene'}));
    resultfig=resultfig(1); % Take the first if there are more.
    set(0,'CurrentFigure',resultfig);

    % get app data
    atlht = getappdata(resultfig,'atlht'); % atlas toolbar
    atlassurfs = getappdata(resultfig,'atlassurfs');
    colorbuttons = getappdata(resultfig,'colorbuttons');

    addht = getappdata(resultfig,'addht'); % object toolbar; see ea_addobj
    idx = zeros(size(addht.Children));
    for i=1:length(addht.Children)
        if strcmp(addht.Children(i).Tag,'roi')
            idx(i) = 1;
        else
            idx(i)=0;
        end
    end
    objects = addht.Children(logical(idx));


    % parse type and toggle
    switch type
        case 'all'
            switch toggle
                case 'on'
                    % atlassurfs
                    idx = [];
                    for i = 1:length(atlassurfs)
                        atlasTag = lower(atlassurfs{i}.Tag);
                        atlasName = regexprep(atlasTag, '_(left|right|midline|mixed)$', '');
                        if strcmp(toggle, 'on') ... % Turn on all atlases
                                || ismember(atlasTag, toggle) ... % Match exactly
                                || ismember(atlasName, toggle) % Match name regardless of side
                            idx = [idx, i];
                            atlassurfs{i}.Visible = 'on';
                            colorbuttons(i).State = 'on';
                        end
                    end
                    atlassurfs = atlassurfs(idx);
                    % objects
                    for i = 1:length(objects)
                        objects(i).State='on';
                    end
                case 'off'
                    % atlassurfs
                    for i = 1:length(atlassurfs)
                        atlassurfs{i}.Visible = 'off';
                        colorbuttons(i).State = 'off';
                    end
                    % objects
                    for i = 1:length(objects)
                        objects(i).State='off';
                    end
            end

        case {'atlas','atlases','atlassurfs'}
            switch toggle
                case 'on'
                    % turn atlassurfs on
                    idx = [];
                    for i = 1:length(atlassurfs)
                        atlasTag = lower(atlassurfs{i}.Tag);
                        atlasName = regexprep(atlasTag, '_(left|right|midline|mixed)$', '');
                        if strcmp(toggle, 'on') ... % Turn on all atlases
                                || ismember(atlasTag, toggle) ... % Match exactly
                                || ismember(atlasName, toggle) % Match name regardless of side
                            idx = [idx, i];
                            atlassurfs{i}.Visible = 'on';
                            colorbuttons(i).State = 'on';
                        end
                    end
                    atlassurfs = atlassurfs(idx);
                    % turn atlassurfs off
                case 'off'
                    for i = 1:length(atlassurfs)
                        atlassurfs{i}.Visible = 'off';
                        colorbuttons(i).State = 'off';
                    end
            end

        case {'object','objects','roi','rois'}
            switch toggle
                case 'on'
                    for i = 1:length(objects)
                        objects(i).State='on';
                    end
                case 'off'
                    for i = 1:length(objects)
                        objects(i).State='off';
                    end
            end


    end
catch

    % LEGACY
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
end