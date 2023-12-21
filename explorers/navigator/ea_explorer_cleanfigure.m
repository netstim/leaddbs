function ea_explorer_cleanfigure(obj)
%% remove existing streamlines
if isempty(obj.drawnstreamlines) % check if prior object has been stored
    if isprop(obj.resultfig,['dt_',obj.ID])
        obj.drawnstreamlines=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.    
    end
end
for i=1:numel(obj.drawnstreamlines)
    try
        delete(obj.drawnstreamlines{i});
    end
end
%% remove existing sweetspots
if isempty(obj.drawnsweetspots) % check if prior object has been stored
    if isprop(obj.resultfig,['dt_',obj.ID])
        obj.drawnsweetspots=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
    end
end
if isfield(obj.drawnsweetspots,'pos')
    for i=1:numel(obj.drawnsweetspots.pos)     
        if ~isempty(obj.drawnsweetspots.pos{i})
            delete(obj.drawnsweetspots.pos{i}.toggleH)
            delete(obj.drawnsweetspots.pos{i})
            obj.drawnsweetspots.pos{i}=[]; 
        end
    end
end
if isfield(obj.drawnsweetspots,'neg')
    for i=1:numel(obj.drawnsweetspots.neg)
        if ~isempty(obj.drawnsweetspots.neg{i})
            delete(obj.drawnsweetspots.neg{i}.toggleH)
            delete(obj.drawnsweetspots.neg{i})
            obj.drawnsweetspots.neg{i}=[];    
        end
    end
end
%% clear remaining buttons
addht=getappdata(obj.resultfig, 'addht');
if ~isempty(addht) && ~isempty(addht.Children)
    delete(addht.Children)    
end
end