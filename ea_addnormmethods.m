function ea_addnormmethods(handles,options,handlestring)

if ~exist('handlestring','var')
    handlestring='normmethod';
end

% add normalization methods to menu
cnt=1;

earoot=ea_getearoot;
ndir=dir([earoot,'ea_normalize_*.m']);

for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);    
        [thisndc,compat]=eval([methodf,'(','''prompt''',')']);
        if compat
            ndc{cnt}=thisndc;
            normmethod{cnt}=methodf;
            try % default function could be something nonexistent
                if strcmp(ndc{cnt},eval([options.prefs.normalize.default,'(','''prompt''',')']))
                    defentry=cnt;
                end
            end
                cnt=cnt+1;
        end
end

try
    setappdata(handles.leadfigure,'normmethod',normmethod);
    set(handles.(handlestring),'String',ndc);
catch
    if isempty(which('spm'))
    warning('It seems that SPM is not installed.');
    end
end

try % set selection of normmethod to default entry (specified in ea_prefs).
    if defentry<=length(get(handles.(handlestring),'String'))
        set(handles.(handlestring),'Value',defentry);
    end
end

clear defentry

ea_switchnormmethod(handles,handlestring);
