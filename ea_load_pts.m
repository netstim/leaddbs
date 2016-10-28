function ea_load_pts(handles,uipatdir,patsub)

if ~exist('patsub','var')
    patsub='patients';
end

if length(uipatdir)>1
    set(handles.patdir_choosebox,'String',['Multiple (',num2str(length(uipatdir)),')']);
    set(handles.patdir_choosebox,'TooltipString',ea_strjoin(uipatdir,', '));
else
    set(handles.patdir_choosebox,'String',uipatdir{1});
    set(handles.patdir_choosebox,'TooltipString',uipatdir{1});
end

% store patient directories in figure


setappdata(handles.leadfigure,'uipatdir',uipatdir);
try
ea_switchctmr(handles);
end

ea_getui(handles); % update ui from patient
ea_storeui(handles); % save in pt folder
ea_addrecentpatient(handles,uipatdir,['Recent ',patsub,':'],patsub);



% add VATs to seeds for connectome mapper case
if isfield(handles,'seeddefpopup')
    
    for pt=1:length(uipatdir)
    direc=[uipatdir{pt},filesep];
    stims=ea_dir2cell(dir([direc,'stimulations']));
    if ~exist('remstims','var')
        remstims=stims;
    else
        todel=[];
        for d=1:length(remstims)
            if ~ismember(remstims{d},stims)
               todel(end+1)=d; 
            end
        end
        remstims(todel)=[];
    end
    end
    remstims=ea_prependvat(remstims);
    set(handles.seeddefpopup,'String',[{'Manually choose seeds'},remstims]);
    
    
    % update cons
    if ~strcmp(get(handles.patdir_choosebox,'String'),'Choose Patient Directory')
        directory=get(handles.patdir_choosebox,'String');
        [~,ptname]=fileparts(directory);
        selectedparc='nan';
        options.prefs=ea_prefs(ptname);
        [mdl,sf]=ea_genmodlist([directory,filesep],selectedparc,options);
        ea_updatemodpopups(mdl,sf,handles);
    end
    
    
end

function ea_updatemodpopups(mdl,sf,handles)

set(handles.fiberspopup,'String',mdl(sf==1));
set(handles.fmripopup,'String',mdl(sf==2));
if isempty(get(handles.fiberspopup,'String'))
    set(handles.fiberspopup,'String','No structural connectome found.');
end
if isempty(get(handles.fmripopup,'String'))
    set(handles.fmripopup,'String','No functional connectome found.');
end



function remstims=ea_prependvat(remstims)
for rs=1:length(remstims)
    remstims{rs}=['Use VATs: ',remstims{rs}];
end