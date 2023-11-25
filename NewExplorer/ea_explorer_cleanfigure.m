function ea_explorer_cleanfigure(obj)
%% remove existing streamlines
if isempty(obj.drawnstreamlines) % check if prior object has been stored
    obj.drawnstreamlines=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
end
for i=1:numel(obj.drawnstreamlines)
    try
        delete(obj.drawnstreamlines{i});
    end
end
%% remove existing sweetspots
if isempty(obj.drawnsweetspots) % check if prior object has been stored
    obj.drawnsweetspots=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
end
if isfield(obj.drawnsweetspots,'pos')
    for i=1:numel(obj.drawnsweetspots.pos)
        try
            delete(obj.drawnsweetspots.pos{i});
        end
    end
end
if isfield(obj.drawnsweetspots,'neg')
    for i=1:numel(obj.drawnsweetspots.neg)
        try
            delete(obj.drawnsweetspots.neg{i});
        end
    end
end
end