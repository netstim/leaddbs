function ea_addnormmethods(handles,options,mstr)

if ~exist('mstr','var')
    mstr=''; % no macaque mode
end

% add normalization methods to menu
cnt=1;

earoot=ea_getearoot;
ndir=dir([earoot,mstr,'ea_normalize_*.m']);

for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);
    try
        [thisndc,compat]=eval([methodf,'(','''prompt''',')']);
        if compat
            ndc{cnt}=thisndc;
            normmethod{cnt}=methodf;
            if strcmp(ndc{cnt},eval([options.prefs.normalize.default,'(','''prompt''',')']))
                defentry=cnt;
            end
            cnt=cnt+1;
        end
    end
end



try
    setappdata(handles.leadfigure,'normmethod',normmethod);
    set(handles.normmethod,'String',ndc);
catch
    if isempty(which('spm'))
    warning('It seems that SPM is not installed.');
    end
end

try % set selection of normmethod to default entry (specified in ea_prefs).
    if defentry<=length(get(handles.normmethod,'String'))
        set(handles.normmethod,'Value',defentry);
    end
end

clear defentry

ea_switchnormmethod(handles);