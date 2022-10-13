function uipatdir = ea_getdataset(options,app)

uipatdir = ea_uigetdir(ea_startpath, 'Please choose dataset folder...');

if isempty(uipatdir)
    return
end

if exist('handles','var')
    ea_load_pts_dataset(app,uipatdir);

    if isfield(app,'atlassetpopup') % not present in connectome mapper
        atlasset=get(app.atlassetpopup,'String');
        atlasset=atlasset{get(app.atlassetpopup,'Value')};

        ea_listatlassets(options,app,get(app.vizspacepopup,'Value'),atlasset);
    end
end
