function uipatdir = ea_getpatients(options,handles)

uipatdir = ea_uigetdir(ea_startpath, 'Please choose patient folder(s)...');

if isempty(uipatdir)
    return
end

if exist('handles','var')
    ea_load_pts(handles,uipatdir);

    if isfield(handles,'atlassetpopup') % not present in connectome mapper
        atlasset = handles.atlassetpopup.String;
        try
            atlasset = atlasset{handles.atlassetpopup.Value};
        catch
            prefs = ea_prefs;
            atlasset =  prefs.machine.defaultatlas;
        end

        ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);
    end
end
