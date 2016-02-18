function ea_cvshowvatdmri(resultfig,directory,handles,selectedparc,options)

% get everything we need to know from handles:
if isstruct(handles) % called from GUI
    vatmodality=get(handles.vatmodality,'String'); vatmodality=vatmodality{get(handles.vatmodality,'Value')};
    vs=get(handles.vatseed,'String'); % dont need this below
    vss=get(handles.vatseed,'Value'); % dont need this below
    vsname=vs{vss};
    [usevat,dimensionality,~,sides]=ea_checkvatselection(handles);
    thresh=get(handles.vatthresh,'String');
elseif iscell(handles) % called from lead_group
    vatmodality=handles{1};
    vsname=handles{2};
    thresh='auto';
    usevat={'right','left'};
    dimensionality=2; % how many ROI.
    sides=[1,2];
end
    


        % fibers filename
        switch vatmodality
            case 'Patient-specific fiber tracts'
                fibersfile=[directory,options.prefs.FTR_normalized];
            otherwise
                fibersfile=[options.earoot,'fibers',filesep,vatmodality,'.mat'];
        end
        
        % seed filename

        seedfile={};
        for v=1:dimensionality
            
            seedfile{v}=[directory,'stimulations',filesep,vsname,filesep,'vat_',usevat{v},'.nii'];
        end
        for side=sides
try
            load([directory,'stimulations',filesep,vsname,filesep,'stimparameters_',usevat{side},'.mat']);
catch
    keyboard
end
        end
        
        targetsfile=[options.earoot,'templates',filesep,'labeling',filesep,selectedparc,'.nii'];
        
        options.writeoutstats=1;
        options.writeoutpm=0;
        
        
        [changedstates,ret]=ea_checkfschanges(resultfig,fibersfile,seedfile,targetsfile,thresh,'vat');
        
        if ~ret % something has changed since last time.
            ea_deletePL(resultfig,'PL','vat');
            if dimensionality % one of the vat checkboxes is active
                [~,thresh]=ea_cvshowfiberconnectivities(resultfig,fibersfile,seedfile,targetsfile,thresh,sides,options,S,changedstates,'vat'); % 'vat' only used for storage of changes.
                if isstruct(handles)
                    set(handles.vatthreshis,'String',num2str(thresh));
                end
            end
            
        end