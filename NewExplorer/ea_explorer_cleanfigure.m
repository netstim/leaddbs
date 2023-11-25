function ea_explorer_cleanfigure(obj)
%% remove existing streamlines
if isempty(obj.drawnstreamlines) % check if prior object has been stored
    obj.drawnstreamlines=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
end
for i=1:numel(obj.drawnstreamlines)
    % try since could run into error when reopening from scratch.
    try
        delete(obj.drawnstreamlines{i});
    end
end
%% remove existing sweetspots
if isempty(obj.drawnsweetspots) % check if prior object has been stored
    obj.drawnsweetspots=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
end
for i=1:numel(obj.drawnsweetspots)
    % try since could run into error when reopening from scratch.
    try
        delete(obj.drawnsweetspots{i});
    end
end
end