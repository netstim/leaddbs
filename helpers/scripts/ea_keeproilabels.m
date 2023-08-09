function ea_keeproilabels(varargin)
%
% Helper function to keep roi labels or turn them off
%
%
% Example:
%   Default to all objects/atlassurfs in scene
%
%   ea_keeproilabels('on'); 
%   ea_keeproilabels('off');
%
%   Note: this syntax has a similar function to ea_keepatlaslabels 
%       and will toggle all objects/atlassurfs in scene on/off.
%       Use ea_keepatlaslabels to toggle type of atlas labels on/off 
%       e.g. ea_keepatlaslabels('STN')
%
% Example 2: include type of object either before or after 'on'/'off'
%    type of object:
%           1. atlas, atlasees, or atlassurfs
%           2. object, objects, roi, rois
%
%   ea_keeproilabels('off','rois');
%   ea_keeproilabels('off','atlases');
%   ea_keeproilabels('rois','on');
%   ea_keeproilabels('on','atlas');
% _________________________________________________________________________
% Copyright (C) 2023 Brigham and Women's Hospital
%
% Ari Kappel

% Parse inputs
% Default: Show all if no input parameters provided
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
atlassurfs = getappdata(resultfig,'atlassurfs');
colorbuttons = getappdata(resultfig,'colorbuttons');

addht = getappdata(resultfig,'addht');
rois = findobj(addht.Children, 'Type', 'uitoggletool', 'UserData', 'roi');

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
                % rois
                for i = 1:length(rois)
                    rois(i).State='on';
                end
            case 'off'
                % atlassurfs
                for i = 1:length(atlassurfs)
                    atlassurfs{i}.Visible = 'off';
                    colorbuttons(i).State = 'off';
                end
                % rois
                for i = 1:length(rois)
                    rois(i).State='off';
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
                for i = 1:length(rois)
                    rois(i).State='on';
                end
            case 'off'
                for i = 1:length(rois)
                    rois(i).State='off';
                end
        end
end