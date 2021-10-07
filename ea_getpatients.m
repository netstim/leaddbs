function uipatdir = ea_getpatients(options,handles)

uipatdir = ea_uigetdir(ea_startpath, 'Please choose patient folder(s)...');

if isempty(uipatdir)
    return
end

if exist('handles','var')
    ea_load_pts(handles,uipatdir);

    if isfield(handles,'atlassetpopup') % not present in connectome mapper
        atlasset=get(handles.atlassetpopup,'String');
        atlasset=atlasset{get(handles.atlassetpopup,'Value')};

        ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);
    end
end
