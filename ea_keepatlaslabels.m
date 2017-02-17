function [atlassurfs] = ea_keepatlaslabels(varargin)
%
% Small function to keep atlas labels given string
%
% Example:
%
%   ea_keepatlaslabels('STN','RN','GPi','GPe')
%   atlassurfs = ea_keepatlaslabels('STN','RN','GP')
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel


H = findall(0,'type','figure');
resultfig = H(~cellfun(@isempty,strfind({H(:).Name},{'Electrode-Scene'})));
atlassurfs = getappdata(resultfig,'atlassurfs');

set(0,'CurrentFigure',resultfig)
idx=zeros(length(atlassurfs),1);
for i = 1:length(varargin)
    idx = idx+~cellfun(@isempty,strfind(get(atlassurfs(:),'Tag'),varargin{i}));
end
    set(atlassurfs(idx==0),'Visible','off')
    set(atlassurfs(idx>0),'Visible','on')
    
atlassurfs = atlassurfs(idx>0);


% view(0,90) %Z-axis view(180,-90)
% view(90,0) %X-axis
% view(180,0) %Y-axis
