function ea_generatemultiplevats(min_amplitude,max_amplitude,stepsize,domonopolar,M,patselect,selected_connectome,result_fig)

% This function is a modified version of the "calculatebutton_Callback" from 
% lead_group.m... I did not clean it up so there is still a lot of stuff in
% here from the original script which will be unnecessary. The basic
% structure of the script is:

% 1) Load and define options
% 2) Open popup window to choose between retrieving data from excel sheet
% or creating pre-defined monopolar stimulations across all contacts
% 3) Iteration through selected patients
% 4) Iterate through selected stimulation settings

% hObject    handle to calculatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Load and define options
%ea_refresh_lg(handles);
options=ea_setopts_local;
va = 0;             % Constant current stimulation
% va = 1;             % Constant voltage stimulation
options.native = 0;             %% Not sure if native = 1 works...
%stimname=ea_detstimname(options);

options.groupmode = 1;
options.groupid = M.guid;

% determine if fMRI or dMRI - Not sure what happens here ;)
%mods=get(handles.fiberspopup,'String');
%mod=mods{get(handles.fiberspopup,'Value')};
%switch mod
%    case {'Patient''s fiber tracts', 'Patient''s fMRI time courses'}
%        fibersfile=mod;
%    case 'Do not calculate connectivity stats'
%    otherwise % load fibertracts once and for all subs here.
        
%end
mod = selected_connectome;
[fibersfile.fibers,fibersfile.fibersidx]=ea_loadfibertracts([ea_getconnectomebase('dmri'),mod,filesep,'data.mat']);



%% Open popup window


% setappdata(multvt,'selection',selection);
% setappdata(multvt,'list',M.patient.list);
% setappdata(multvt,'amps',[]);
% setappdata(multvt,'multset','none');
% setappdata(multvt,'stimsets',[]);
% uicontrol('Parent',multvt,'Style','pushbutton','String','Simulate multiple monopolar stimulation settings','Position',[10 55 430 35],'FontSize',14,'Callback',@ea_simulatestims);
% uicontrol('Parent',multvt,'Style','pushbutton','String','Extract settings from excel sheet','Position',[10 10 430 35],'FontSize',14,'Callback',@ea_extractstimsfromxls);
% 
% set(multvt, 'Visible', 'on');
% waitfor(multvt,'UserData',1)    % Waiting to continue script after callack was executed

selection = linspace(1,length(patselect),length(patselect)); %selection of patients
amps.min = min_amplitude;
amps.max = max_amplitude;
amps.stepsize = stepsize;
if domonopolar
    multset = 'mon';
end


%% Start iterating through patients
for pt=selection
    
    %% set pt specific options
    % Again a lot of lines which I did not touch (I think you can ignore
    % this section
    
    % own fileparts to support windows/mac/linux slashes even if they come
    % from a different OS.
    if ~strfind(M.patient.list{pt},'/')
        lookfor='\';
    else
        lookfor='/';
    end
    
    slashes=strfind(M.patient.list{pt},lookfor);
    if ~isempty(slashes)
        options.patientname=M.patient.list{pt}(slashes(end)+1:end);
        options.root=M.patient.list{pt}(1:slashes(end));
    else
        options.patientname=M.patient.list{pt};
        options.root='';
    end
    
    fprintf('\nProcessing %s...\n\n', options.patientname);
    try
        options.numcontacts=size(M.elstruct(pt).coords_mm{1},1);
    catch % no localization present or in wrong format.
        ea_error(['Please localize ',options.patientname,' first.']);
    end
    options.elmodel=M.elstruct(pt).elmodel;
    options=ea_resolve_elspec(options);
    options.prefs=ea_prefs(options.patientname);
    options.d3.verbose='off';
    options.d3.elrendering=1;	% hard code to viz electrodes in this setting.
    options.d3.exportBB=0;	% don't export brainbrowser struct by default
    options.d3.colorpointcloud=0;
    options.d3.hlactivecontacts=1;
    options.d3.showactivecontacts =1;
    options.d3.showpassivecontacts=1;
    %options.d3.hlactivecontacts=get(handles.highlightactivecontcheck,'Value');
    %options.d3.showactivecontacts=get(handles.showactivecontcheck,'Value');
    %options.d3.showpassivecontacts=get(handles.showpassivecontcheck,'Value');
    try
        options.d3.isomatrix=M.isomatrix;
    catch
        options.d3.isomatrix={};
    end
    
    options.d3.isovscloud=M.ui.isovscloudpopup;
    options.d3.showisovolume=M.ui.showisovolumecheck;
    options.d3.exportBB=0;
    options.expstatvat.do=0;
    try
        options.expstatvat.vars=M.clinical.vars(M.ui.clinicallist);
        options.expstatvat.labels=M.clinical.labels(M.ui.clinicallist);
        options.expstatvat.pt=pt;
    end
    options.expstatvat.dir=M.ui.groupdir;
    processlocal=0;
    
    if M.ui.detached
        processlocal=1;
        mkdir([M.ui.groupdir,options.patientname]);
        options.root=M.ui.groupdir;
        %    options.patientname='tmp';
        try
            ea_stats=M.stats(pt).ea_stats;
        catch
            ea_stats=struct;
        end
        reco.mni.coords_mm=M.elstruct(pt).coords_mm;
        reco.mni.trajectory=M.elstruct(pt).trajectory;
        reco.mni.markers=M.elstruct(pt).markers;
        reco.props.elmodel=M.elstruct(pt).elmodel;
        reco.props.manually_corrected=1;
        save([M.ui.groupdir,options.patientname,filesep,'ea_stats'],'ea_stats');
        save([M.ui.groupdir,options.patientname,filesep,'ea_reconstruction'],'reco');
    end
    
    if ~exist(options.root,'file') % data is not there. Act as if detached. Process in tmp-dir.
        processlocal=1;
        warning('on');
        warning('Data has been detached from group-directory. Will process locally. Please be aware that you might loose this newly-processed data once you re-attach the single-patient data to the analysis!');
        warning('off');
        mkdir([M.ui.groupdir,options.patientname]);
        options.root=M.ui.groupdir;
        % options.patientname='tmp';
        try
            ea_stats=M.stats(pt).ea_stats;
        catch
            ea_stats=struct;
        end
        reco.mni.coords_mm=M.elstruct(pt).coords_mm;
        reco.mni.trajectory=M.elstruct(pt).trajectory;
        reco.mni.markers=M.elstruct(pt).markers;
        reco.props.elmodel=M.elstruct(pt).elmodel;
        reco.props.manually_corrected=1;
        save([M.ui.groupdir,options.patientname,filesep,'ea_stats'],'ea_stats');
        save([M.ui.groupdir,options.patientname,filesep,'ea_reconstruction'],'reco');
    end
    
    %delete([options.root,options.patientname,filesep,'ea_stats.mat']);
    
    options.leadprod = 'group';
    options.patient_list=M.patient.list;
    options.d3.mirrorsides=0;
    resultfig=ea_elvis(options,M.elstruct(pt));
    
    
    % save scene as matlab figure
    
    options.modality=ea_checkctmrpresent(M.patient.list{pt});
    volumespresent=1;
    if options.modality(1) % prefer MR
        options.modality=1;
    else
        if options.modality(2)
            options.modality=2;
        else
            options.modality=1;
            warning(['No MR or CT volumes found in ',M.patient.list{pt},'.']);
            volumespresent=0;
        end
    end
    
    %% Define stimulation settings
    
    if strcmp(multset,'mon')    % True if you chose to calculate multiple monopolar stimulation settings
        
        ncnttmp = options.numcontacts;
        
        allstims{1} = repmat(((1:ncnttmp)-1)',1,size(amps.min:amps.stepsize:amps.max,2));
        allstims{2} = repmat(amps.min:amps.stepsize:amps.max,ncnttmp,1);
        
        runs = 1:numel(allstims{1});
        
    else                        % true if you chose to retrieve stimulation settings from excel sheet
        runs = 1;
    end

    %% Iterate through all stimulations
    if isfield(M,'S')
        for n = runs
            try % Generate M struct according to monopolar settings
                stimtmp = zeros(1,9);
                stimtmp(1,1) = allstims{1,2}(n);
                stimtmp(1,allstims{1,1}(n)+2) = -100;     
                M.S(pt) = jr_generateMfile(stimtmp,stimtmp,M.S(pt),va);
                M.S(pt).label = ['c',num2str(allstims{1}(n),'%02d'),'_a',num2str(allstims{2}(n)*10,'%02d')];
            catch
                ea_error(['Stimulation parameters for ',M.patient.list{pt},' are missing.']);
            end
            
            ea_genvat = eval('@ea_genvat_horn_modified');
            
            
            %% More options
            % Not sure why I this is defined down here and not in the
            % beginning of the script. Also since the patient space
            % calculations were just implemented while I was working on
            % this I needed to do some weird workaround. But since you are
            % calculating everything in template space you should be able
            % to just leave it like this
            options.orignative=options.native; % backup
            options.native=~ea_getprefs('vatsettings.estimateInTemplate'); % see whether VTAs should be directly estimated in template space or not
            options.native = 0;
            if options.native && ~volumespresent
                warning(['You chose to process VTAs in native space but patient-data cannot be found for ',M.patient.list{pt},'. Proceeding with VTA calculation directly in template space.']);
                options.native=0;
            end
            
%             options.prefs.machine.vatsettings.horn_atlasset = options.atlasset;          % Added additionally otherwise DISTAL Atlas will be loaded in headmodel 
            options.prefs.machine.vatsettings_horn_removeElectrode = 1;
            
            % Atlas definitions! I remember I had tried this with the
            % DISTAL atlas as well but gave up at some point because
            % generating the VTAs failed in 10% of the times because of
            % self intersections in the head model. This was a common bug
            % back then not sure it was fixed by now. Just be careful when
            % you select different atlases and check the efields (you
            % should see the outline of the anatomical structures in the
            % efield)
            if ~strcmp(options.atlasset,'Homogeneous')
                options.atlasset = 'Homogeneous';
            else
                options.prefs.machine.vatsettings.horn_cgm = 0.20;
                options.prefs.machine.vatsettings.horn_cwm = 0.20;
                options.prefs.machine.vatsettings.horn_useatlas = 0;
            end
                
            % Define the name of the folder for the nii to be saved in
            if va == 0
                volcur = 'mA';
            elseif va == 1
                volcur = 'V';
            else
                volcur = '??';
            end
            
            if options.prefs.machine.vatsettings_horn_removeElectrode
                elmov = 'elmov';                % this variable defines the "electrode removal" option 
            else
                elmov = 'noelmov';
            end
            
            fname = [elmov,'_', volcur, '_', options.atlasset, '_', num2str(round(options.prefs.machine.vatsettings.horn_cgm*100),'%02d'), '_', num2str(round(options.prefs.machine.vatsettings.horn_cwm*100),'%02d')];
            
            %% Iterate both hemispheres
            for side=1:2
                disp([' ', newline, 'Patient ', M.patient.list{pt}, newline,'Simulating efield: ', fname, ' side ', num2str(side),' | ', M.S(pt).label])
                setappdata(resultfig,'elstruct',M.elstruct(pt));
                setappdata(resultfig,'elspec',options.elspec);
                try
                    feval(ea_genvat,M.elstruct(pt).coords_mm,M.S(pt),side,options,fname,options.prefs.machine.vatsettings.horn_ethresh,result_fig,1,n);
                catch
                    msgbox({['Error while creating VTA of ',M.patient.list{pt}];['Setting: ', fname, ' side ', num2str(side),' | ', M.S(pt).label,'.']});
                end
                %                 stimparams(1,side).volume=volume;
            end
        end
        options.native=options.orignative; % restore
    end
    
end
%% processing done here.

%ea_refresh_lg(handles);


function options=ea_setopts_local

options.earoot=ea_getearoot;
options.verbose=3;
options.sides=1:2; % re-check this later..
options.atlasset = 'Homogeneous';
options.fiberthresh=1;
options.writeoutstats=1;
options.writeoutpm=1;
options.colormap=jet;
options.d3.write=1;
options.d3.prolong_electrode=2;
options.d3.writeatlases=1;
options.macaquemodus=0;
%try
%    options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
%catch % too many entries..
%    set(handles.atlassetpopup,'Value',1);
%    options.atlasset=1;
%end

% options.labelatlas=get(handles.labelpopup,'String');
% try
%     options.labelatlas=options.labelatlas{get(handles.labelpopup,'Value')};
% catch % too many entries..
%     set(handles.labelpopup,'Value',1);
%     options.labelatlas=1;
% end


function multvt = ea_extractstimsfromxls(multvt,~,~)

% rows of xls file specify patients, columns should be structured like in
% this example: leaddbs\helpers\getstimsfromxls_exampl.xlsx
% Tab needs to have the name 'stimsets' and names of patients in column 1
% has to match the folder names of the lead locations
[fname,folder] = uigetfile('*.xls*','Select a file');
[stimsets,names] = xlsread([folder fname], 'stimsets','A3:S95');
selection_old = getappdata(multvt.Parent,'selection');
list = getappdata(multvt.Parent,'list');

ptfound = zeros(length(selection_old),length(names));
for i = 1:length(names)
    ptfound(:,i) = contains(list,names{i});
end

[selection_new,ptselxls] = ind2sub(size(ptfound),find(ptfound));
stimsets = stimsets(ptselxls,:);
ptfound = sum(ptfound,2);

waitfor(msgbox([num2str(sum(double(ptfound),'all')), '/' num2str(length(selection_old)), ' patients were found in excel list and are being processed']));
setappdata(multvt.Parent,'selection',selection_new');
setappdata(multvt.Parent,'multset','xls');
setappdata(multvt.Parent,'stimsets',stimsets);
multvt.Parent.UserData = 1;



function ea_simulatestims(multvt,min_amplitude,max_amplitude,stepsize)
%amps.if = inputdlg({'Min amplitude','Stepsize','Max amplitude'},'Monopolar Review boundaries',[1 35;1 35;1 35],{'5','0.5','5'});
amps.min = min_amplitude;
amps.stepsize = stepsize;
amps.max = max_amplitude;
% amps.min = str2double(amps.if{1});
% amps.stepsize = str2double(amps.if{2});
% amps.max = str2double(amps.if{3});

setappdata(multvt.Parent,'multset','mon');
setappdata(multvt.Parent,'amps',amps);
multvt.Parent.UserData = 1;




