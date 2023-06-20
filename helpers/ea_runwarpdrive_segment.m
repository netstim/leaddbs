function [] = ea_runwarpdrive_segment(~,~,handles)

options = ea_handles2options(handles);
options.uipatdirs = getappdata(handles.leadfigure,'uipatdir');
options.leadprod = 'dbs';
setappdata(handles.leadfigure,'handles',handles);
options.leadfigure=handles.leadfigure;

for i = 1:length(options.uipatdirs)

    options.subjInd=i;

    if isfield(options, 'leadfigure')
        bids = getappdata(options.leadfigure, 'bids');
        subjId = getappdata(options.leadfigure, 'subjId');
        if ~isempty(bids)
            options.bids = bids;
            if ~isempty(subjId{options.subjInd})
                options.subj = bids.getSubj(subjId{options.subjInd}, options.modality);
            end
        end
    end

    options.overwriteapproved = 1;
    ea_runwarpdrive(options, '1');

end


