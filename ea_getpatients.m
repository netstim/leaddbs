function uipatdir=ea_getpatients(options,handles)


p='/'; % default use root
try
    p=pwd; % if possible use pwd instead (could not work if deployed)
end
try % finally use last patient parent dir if set.
    earoot=ea_getearoot;
    load([earoot,'common',filesep,'ea_recentpatients.mat']);
    p=fileparts(fullrpts{1});
end

uipatdir=ea_uigetdir(p,'Please choose patient folder(s)...');

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