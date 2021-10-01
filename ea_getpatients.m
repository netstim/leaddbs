function uipatdir = ea_getpatients(options,handles)

startFolder = ea_gethome; % default home folder

try
    startFolder = pwd; % if possible use pwd instead (could not work if deployed)
end

try % finally use last patient parent dir if set.
    load([ea_getearoot,'common', filesep, 'ea_recentpatients.mat'], 'recentfolders');
    startFolder = fileparts(recentfolders{1});
end

uipatdir = ea_uigetdir(startFolder, 'Please choose patient folder(s)...');

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
