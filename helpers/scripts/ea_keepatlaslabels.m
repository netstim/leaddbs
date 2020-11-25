function [atlassurfs] = ea_keepatlaslabels(varargin)
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

% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel

try isvalid(varagin{1})
    if ~isvalid(varargin{1})
        ea_warning('Figure handle not valid')
    end
catch
    H = findall(0,'type','figure');
    resultfig = H(contains({H(:).Name},{'Electrode-Scene'}));
end

if isempty(varargin) || isempty(varargin{1}) || strcmpi(varargin{1},'on') || ...
        ( length(varargin)==2 && isempty(varargin{2}) )
    varargin{1}='right';
    varargin{2}='left';
end

resultfig=resultfig(1); % take the first if there are more.

atlassurfs = getappdata(resultfig,'atlassurfs');
colorbuttons = getappdata(resultfig,'colorbuttons');
set(0,'CurrentFigure',resultfig)
idx=zeros(length(atlassurfs),1);
for i = 1:length(varargin)
    idx = idx+ismember(get(atlassurfs(:),'Tag'),[varargin{i},'_left']);
    idx = idx+ismember(get(atlassurfs(:),'Tag'),[varargin{i},'_right']);
    idx = idx+ismember(get(atlassurfs(:),'Tag'),[varargin{i},'_midline']);
    idx = idx+ismember(get(atlassurfs(:),'Tag'),[varargin{i},'_mixed']);
end
    set(colorbuttons(idx==0),'State','off')
    set(atlassurfs(idx==0),'Visible','off')
    
    set(colorbuttons(idx>0),'State','on')
    set(atlassurfs(idx>0),'Visible','on')
    
atlassurfs = atlassurfs(idx>0);


% view(0,90) %Z-axis view(180,-90)
% view(90,0) %X-axis
% view(180,0) %Y-axis view(0,0) (COR)
