function ea_init_coregctpopup(handles,options,handlestring)

if ~exist('handlestring','var')
    handlestring='coregctmethod';
end

% add coreg methods to menu
cnt=1;
ndir=dir([ea_getearoot,'ea_coregctmri_*.m']);
for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);
    try
        [thisndc,spmvers]=eval([methodf,'(','''prompt''',')']);
        if ismember(spm('ver'),spmvers)
            cdc{cnt}=thisndc;
            coregctmethod{cnt}=methodf;
            if strcmp(cdc{cnt},eval([options.prefs.ctcoreg.default,'(','''prompt''',')']))
                defentry=cnt;
            end
            cnt=cnt+1;
        end
    catch
        keyboard
    end
end

try
    setappdata(handles.leadfigure,'coregctmethod',coregctmethod);
    set(handles.(handlestring),'String',cdc);
catch
    if isempty(which('spm'))
        ea_error('Please install SPM12 for Lead-DBS to work properly.');
    end
end

try % set selection of ctcoregmethod to default entry (specified in ea_prefs).
    if defentry<=length(get(handles.(handlestring),'String'))
        set(handles.(handlestring),'Value',defentry);
    end
end
