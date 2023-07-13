function ea_refresh_lg(handles)

ea_busyaction('on',handles.leadfigure,'group');

options = getappdata(handles.leadfigure, 'options');

% get model data
disp('Getting model data...');
M=getappdata(handles.leadfigure,'M');

if strcmp(handles.groupdir_choosebox.String,'Choose Dataset Directory') % not set yet.
    ea_busyaction('off',handles.leadfigure,'group');
    return
end

M.root = fullfile(handles.groupdir_choosebox.String,filesep);

if ~isfield(M,'guid') % only done once, legacy support.
    M.guid=datestr(datevec(now), 'yyyymmddHHMMSS' );
end

disp('Refreshing group list...');

% refresh group list
set(handles.grouplist,'String',M.patient.group);
if max(M.ui.listselect) > length(M.patient.group)
    M.ui.listselect = 1;
end
try set(handles.grouplist,'Value',M.ui.listselect); end

disp('Refreshing patient list...');
% refresh patient list
set(handles.patientlist,'String',M.patient.list);
try set(handles.patientlist,'Value',M.ui.listselect); end

disp('Creating isomatrix from regressor list...');
% set isomatrix from variable in clinical list
try
    M.isomatrix=M.clinical.vars(get(handles.clinicallist,'Value'));
    M.isomatrix_name=M.clinical.labels(get(handles.clinicallist,'Value'));
catch
    M.isomatrix={};
    M.isomatrix_name={};
end

disp('Refreshing selections on VI / FC Lists...');

parcellations = get(handles.labelpopup,'String');
if ischar(parcellations)
    parcellations = {parcellations};
end

if ~isfield(M.ui,'labelpopup')
    M.ui.labelpopup = parcellations{get(handles.labelpopup,'Value')};
else
    if isnumeric(M.ui.labelpopup)
        if M.ui.labelpopup>0 && M.ui.labelpopup<=length(parcellations)
            set(handles.labelpopup,'Value',M.ui.labelpopup);
        else % Set to default parcellation in case index out of range
            defaultParc = options.prefs.lg.defaultParcellation;
            set(handles.labelpopup,'Value',find(ismember(parcellations, defaultParc)));
        end
    else
        parcInd = find(ismember(parcellations, M.ui.labelpopup), 1);
        if ~isempty(parcInd)
            set(handles.labelpopup,'Value',parcInd);
        else
            defaultParc = options.prefs.lg.defaultParcellation;
            set(handles.labelpopup,'Value',find(ismember(parcellations, defaultParc)));
        end
    end
end

thisparc = parcellations{get(handles.labelpopup,'Value')};

% dMRI:
if ~isempty(M.patient.list)
    patientDirectory = [M.patient.list{1},filesep];
else
    patientDirectory = [];
end
modlist = ea_genmodlist(patientDirectory,thisparc,options,'dmri');
if ~ismember('Patient''s fiber tracts' ,modlist)
    modlist{end+1}='Patient''s fiber tracts';
end
modlist{end+1}='Do not calculate connectivity stats';
set(handles.fiberspopup,'String',modlist);
if get(handles.fiberspopup,'Value') > length(modlist)
    set(handles.fiberspopup,'Value',1);
end

% update UI
disp('Updating UI...');

disp('Getting stimulation parameters...');
S=getappdata(handles.leadfigure,'S');
S=ea_checkS(M,S,options,handles);
if ~isempty(S)
    set(handles.setstimparamsbutton,'BackgroundColor',[0.1;0.8;0.1]);
    M.S=S;
    for sp=1:length(M.S) % make sure stimlabel is gs_guid.
        M.S(sp).label=['gs_',M.guid];
    end
    % M.S=ea_activecontacts(M.S);
    M.vatmodel=getappdata(handles.leadfigure,'vatmodel');
else
    set(handles.setstimparamsbutton,'BackgroundColor',[0.93,0.93,0.93]);
end

% check if groups are okay
if isfield(M,'groups')
    if ~isequal((unique(M.patient.group)),M.groups.group)
        % reassign groups and colors
        M.groups.group=unique(M.patient.group);
        C=ea_color_wes('all');
        C=rgb2hsv(C);
        C(:,2)=C(:,2)./2;
        C=hsv2rgb(C);
        M.groups.color=C(M.groups.group,:);
        M.groups.colorschosen=1;
    end
    % make choosecolors button green if chosen.
    if isfield(M.groups,'colorschosen')
        set(handles.choosegroupcolors,'BackgroundColor',[0.1;0.8;0.1]);
    else
        set(handles.choosegroupcolors,'BackgroundColor',[0.93,0.93,0.93]);
    end
end

% update checkboxes:
try set(handles.showactivecontcheck,'Value',M.ui.showactivecontcheck); end
try set(handles.showpassivecontcheck,'Value',M.ui.showpassivecontcheck); end
try set(handles.highlightactivecontcheck,'Value',M.ui.hlactivecontcheck); end
try set(handles.showisovolumecheck,'Value',M.ui.showisovolumecheck); end
try set(handles.statvatcheck,'Value',M.ui.statvat); end

% update atlas selectboxes
atlasset = get(handles.atlassetpopup,'String');
if ~isfield(M.ui,'atlassetpopup')
    M.ui.atlassetpopup = atlasset{get(handles.atlassetpopup,'Value')};
else
    if isnumeric(M.ui.atlassetpopup) % Compatible with old lead group file
        if M.ui.atlassetpopup>0 && M.ui.atlassetpopup<=length(atlasset)
            set(handles.atlassetpopup,'Value',M.ui.atlassetpopup);
        else % Set to default atlas in case index out of range
            defaultAtlas = options.prefs.machine.defaultatlas;
            set(handles.atlassetpopup,'Value',find(ismember(atlasset, defaultAtlas)));
        end
    else % New lead group file in which atlassetpopup is the name of the atlas
        atlasInd = find(ismember(atlasset, M.ui.atlassetpopup), 1);
        if ~isempty(atlasInd)
            set(handles.atlassetpopup,'Value',atlasInd);
        else
            defaultAtlas = options.prefs.machine.defaultatlas;
            set(handles.atlassetpopup,'Value',find(ismember(atlasset, defaultAtlas)));
        end
    end
end

connectomes = get(handles.fiberspopup,'String');
defaultConnectomeIdx = length(connectomes);
if ~(ischar(connectomes) && strcmp(connectomes, 'Fibers'))
    if ~isfield(M.ui, 'connectomename')
        set(handles.fiberspopup,'Value',defaultConnectomeIdx);
        M.ui.connectomename = connectomes{defaultConnectomeIdx};
    else
        connectomeIdx = find(ismember(connectomes, M.ui.connectomename),1);
        if ~isempty(connectomeIdx)
            set(handles.fiberspopup,'Value',connectomeIdx);
        else
            set(handles.fiberspopup,'Value',defaultConnectomeIdx);
        end
    end
end

try set(handles.normregpopup,'Value',M.ui.normregpopup); end

t=datetime('now');
t.Format='uuuMMddHHmmss';
t=str2double(char(t));
if ~isfield(M.ui,'lastupdated') || t-M.ui.lastupdated>0 % 0 mins time limit
    % patient specific part:
    if ~isempty(M.patient.list)
        disp('Loading localizations...');
        for pt=1:length(M.patient.list)
            % set stimparams based on values provided by user
            for side=1:2
                if M.ui.labelpopup>length(get(handles.labelpopup,'String'))
                    M.ui.labelpopup=length(get(handles.labelpopup,'String'));
                end
            end
            
            if isfield(M,'stimparams') % deprecated.
                M=rmfield(M,'stimparams');
            end
            
            % load localization
            [~, patientname] = fileparts(M.patient.list{pt});
            
            M.elstruct(pt).group=M.patient.group(pt);
            if ~isfield(M,'groups')
                M.groups.color=[0.7,0.7,0.7];
                M.groups.group=ones(length(M.patient.list),1);
            end
            
            M.elstruct(pt).groupcolors=M.groups.color;
            M.elstruct(pt).groups=M.groups.group;
            
            options.sides=1:2;
            options.native=0;
            try
                [options.root, options.patientname] = fileparts(M.patient.list{pt});
                options.root = [options.root, filesep];
                
                options.subj.recon.recon = fullfile(options.root, options.patientname, 'reconstruction', [options.patientname, '_desc-reconstruction.mat']);
                options = ea_resolve_elspec(options);

                if isfile(options.subj.recon.recon)
                    [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(options);
                    
                    if M.ui.elmodelselect==1 % use patient specific elmodel
                        if exist('elmodel','var')
                            M.elstruct(pt).elmodel=elmodel;
                        else
                            M.elstruct(pt).elmodel='Medtronic 3389'; % use default for older reconstructions that did not store elmodel.
                        end
                    else
                        elmodels = [{'Patient specified'}; ea_resolve_elspec];
                        M.elstruct(pt).elmodel = elmodels{M.ui.elmodelselect};
                    end
                    
                    % make sure coords_mm is congruent to coded electrode model
                    poptions=options;
                    poptions.native=0;
                    poptions.elmodel=M.elstruct(pt).elmodel;
                    poptions=ea_resolve_elspec(poptions);
                    [coords_mm,trajectory,markers]=ea_resolvecoords(markers,poptions,0);
                    
                    M.elstruct(pt).coords_mm=coords_mm;
                    M.elstruct(pt).coords_acpc=coords_acpc;
                    M.elstruct(pt).trajectory=trajectory;
                    M.elstruct(pt).name = patientname;
                    
                    if ~exist('markers','var') % backward compatibility to old recon format
                        for side=1:2
                            try
                                %if side is present, process the old recon format
                                markers(side).head=coords_mm{side}(1,:);
                                markers(side).tail=coords_mm{side}(4,:);
                                [xunitv, yunitv] = ea_calcxy(markers(side).head, markers(side).tail);
                                markers(side).x = markers(side).head +  xunitv*(options.elspec.lead_diameter/2);
                                markers(side).y = markers(side).head + yunitv*(options.elspec.lead_diameter/2);
                            catch
                            end
                        end
                    end
                    M.elstruct(pt).markers=markers;
                end
            catch
                if pt>1 % first patient has worked but some other patient seems not to have worked.
                    try
                        M.elstruct(1).coords_mm; % probe if error happens in pt. 1 ? if not show warning
                        warning(['No reconstruction present for ',patientname,'. Please check.']);
                    end
                end
            end
        end
        
        %uniform the data (by checking the missing sides and filling them)
        num_sides=length(options.sides);%minimum number of sides is 2 (R and L); (Hardcorded for now)
        for pt=1:length(M.patient.list)
            if length(M.elstruct(pt).coords_mm)>num_sides
                num_sides=M.elstruct(pt).coords_mm;
            end
        end
        for pt=1:length(M.patient.list)
            for check_side=1:num_sides %options.sides
                if ea_arenopoints4side(M.elstruct(pt).coords_mm, check_side)
                    %force to have empty values if side is not present
                    M.elstruct(pt).coords_mm{check_side}=[];
                    if isfield(M.elstruct(pt),'coords_acpc') && iscell(M.elstruct(pt).coords_acpc)
                        if ea_arenopoints4side(M.elstruct(pt).coords_acpc, check_side)
                            M.elstruct(pt).coords_acpc{check_side}=[];
                        end
                    end
                    M.elstruct(pt).trajectory{check_side}=[];
                    
                    %this will create the missing structure
                    M.elstruct(pt).markers(check_side).head=[];
                    M.elstruct(pt).markers(check_side).tail=[];
                    M.elstruct(pt).markers(check_side).x=[];
                    M.elstruct(pt).markers(check_side).y=[];
                end
            end
        end
        
        % load stats for group
        disp('Loading stats for group...');
        for pt=1:length(M.patient.list)
            % (re-)load stats
            try
                [~, patientname] = fileparts(M.patient.list{pt});
                statsFile = [M.patient.list{pt}, filesep, patientname, '_desc-stats.mat'];
                load(statsFile, 'ea_stats');
                ea_stats=ea_rmssstimulations(ea_stats,M); % only preserve stimulations with label 'gs_groupid'.
                M.stats(pt).ea_stats=ea_stats;
                if isfield(M.stats(pt).ea_stats.atlases,'rebuild') % old stats format with complete atlas table - delete, will lead to large M file
                    M.stats(pt).ea_stats=rmfield(M.stats(pt).ea_stats,'atlases');
                    M.stats(pt).ea_stats.atlases.names=ea_stats.atlases.names;
                    M.stats(pt).ea_stats.atlases.types=ea_stats.atlases.types;
                    
                    % also correct single subject file:
                    load(statsFile, 'ea_stats');
                    ea_stats.atlases=M.stats(pt).ea_stats.atlases;
                    save(statsFile, 'ea_stats', '-v7.3');
                end
            catch ME
                ea_cprintf('CmdWinWarnings', '%s\n', ME.message);
            end
            
            if ~isfield(M,'stats')
                % if no stats  present yet, return.
                setappdata(handles.leadfigure,'M',M);
                set(handles.leadfigure,'name','Lead Group Analysis');
                set(handles.calculatebutton, 'BackgroundColor', [0.1;0.8;0.1]);
                set(handles.explorestats, 'Enable', 'off');
                set(handles.exportstats, 'Enable', 'off');
                break
            else
                set(handles.calculatebutton, 'BackgroundColor', [0.93,0.93,0.93]);
                set(handles.explorestats, 'Enable', 'on');
                set(handles.exportstats, 'Enable', 'on');
            end
            
            priorvilist=M.vilist;
            try % try using stats from patient folder.
                M.vilist=ea_stats.atlases.names;
            catch
                try % try using stats from M-file.
                    M.vilist=M.stats(pt).ea_stats.atlases.names;
                catch
                    M.vilist={};
                end
            end
            
            %disp('Comparing stats with prior atlas intersection list...');
            % check and compare with prior atlas intersection list.
            if ~isempty(priorvilist) && ~isequal(priorvilist,M.vilist)
                warning('off', 'backtrace');
                [~, ptname] = fileparts(M.patient.list{pt});
                warning('%s: inhomogeneous stats. Please re-run group analysis (Calculate Stats).', ptname);
                warning('on', 'backtrace');
            end
            
            priorfclist=M.fclist;
            try % try using stats from patient folder.
                M.fclist=ea_stats.stimulation(1).ft(1).labels{1};
                fcdone=1;
            catch
                try % try using stats from M-file.
                    M.fclist=M.stats(pt).ea_stats.stimulation(1).ft(1).labels{1};
                    fcdone=1;
                catch
                    M.fclist={};
                    fcdone=0;
                end
            end
            
            % check and compare with prior fibertracking list.
            if fcdone
                if ~isempty(priorfclist) && ~isequal(priorfclist,M.fclist)
                    warning('Trying to analyse inhomogeneous patient group. Please re-run single subject lead analysis with patients using always the same labeling atlas.');
                end
            end
        end
        
        
        % load clinical data for group
        disp('Loading clinical data for group...');
        for pt=1:length(M.patient.list)
            
            if exist(fullfile(M.patient.list{pt},'clinical','clinical_scores.mat'),'file')
                ptscores=load(fullfile(M.patient.list{pt},'clinical','clinical_scores.mat'));
                if ~exist('ea_scores','var') % set up automatic assignment struct
                    ea_scores=load(fullfile(ea_getearoot,'clinical','ea_scores.mat'));
                end
                entries=fieldnames(ptscores);
                for entry=1:length(entries) % iterate scores available in patient
                    
                    switch entries{entry}
                        case {'Motor_UPDRS','Motor_MDSUPDRS'} % group these two together
                            
                            scorename='mUPDRS';
                        otherwise
                            
                            scorename=entries{entry};
                    end
                    if exist('clindata','var') && isfield(clindata,scorename)
                    else
                        clindata.(scorename).baseline=[];
                        clindata.(scorename).postop=[];
                        clindata.(scorename).factors={};
                        clindata.(scorename).score={};
                        clindata.(scorename).somatotopies=struct;
                        clindata.(scorename).somatotopynames={};
                    end
                    
                    [ispresent,ix]=ismember(entries{entry},fieldnames(ea_scores));
                    if ispresent % does lead-dbs know what to do with a score named like this?
                        if isfield(ptscores.(entries{entry}),ea_scores.(entries{entry}).default_baseline)
                            baseline=ptscores.(entries{entry}).(ea_scores.(entries{entry}).default_baseline).score;
                            success=1;
                        else
                            success=0;
                            continue
                        end
                        if isfield(ptscores.(entries{entry}),M.guid)
                            postop=ptscores.(entries{entry}).(M.guid).score;
                            success=1;
                        else
                            success=0;
                            continue
                        end
                        if success % both baseline and postop available
                            clindata.(scorename).baseline(pt,:)=table2array(baseline);
                            clindata.(scorename).postop(pt,:)=table2array(postop);
                            clindata.(scorename).score{pt,1}=entries{entry};
                            clindata.(scorename).factornames=ea_scores.(entries{entry}).factornames;
                            clindata.(scorename).factors.(entries{entry})=ea_scores.(entries{entry}).factors;
                            clindata.(scorename).somatotopies.(entries{entry})=ea_scores.(entries{entry}).somatotopies;
                            clindata.(scorename).somatotopynames=ea_scores.(entries{entry}).somatotopynames;
                        end
                    end
                    
                end
            end
        end
        if exist('clindata','var')
            fns=fieldnames(clindata);
            for fn=1:length(fns)
                scorename=fns{fn};
                clindata.(scorename).scores=unique(clindata.(scorename).score);
                
                modes={'absolute','percent','cleaned'};
                
                % now auto-query improvements for all factors:
                for s=1:length(clindata.(scorename).somatotopynames)
                    for f=1:length(clindata.(scorename).factornames)
                        for mode=1:length(modes)
                            I=ea_getimprovs_fctr_smtp(clindata.(scorename),modes{mode},1:length(M.patient.list),clindata.(scorename).factornames{f},clindata.(scorename).somatotopynames{s});
                            Ilabel=[scorename,'_',clindata.(scorename).somatotopynames{s},'_',clindata.(scorename).factornames{f},'_',modes{mode}];
                            
                            [is,ix]=ismember(Ilabel,M.clinical.labels);
                            if ~is
                                M.clinical.labels{end+1}=Ilabel;
                                M.clinical.vars{end+1}=I;
                            else % update
                                M.clinical.labels{ix}=Ilabel;
                                M.clinical.vars{ix}=I;
                            end
                        end
                    end
                end
            end
        end
        
        
        % sync stimulation parameters for group
        disp('Syncing stimulation parameters for group...');
        for pt=1:length(M.patient.list)
            stimParamFile = fullfile(M.patient.list{pt}, 'stimulations', ea_nt(0), ['gs_', M.guid], [options.patientname, '_desc-stimparameters.mat']);
            if isfile(stimParamFile)
                ptS=load(stimParamFile);

                try % could fail if M.S(pt) is not defined.
                    if ~all([isequal(ptS.S.Rs1,M.S(pt).Rs1),...
                            isequal(ptS.S.Rs2,M.S(pt).Rs2),...
                            isequal(ptS.S.Rs3,M.S(pt).Rs3),...
                            isequal(ptS.S.Rs4,M.S(pt).Rs4),...
                            isequal(ptS.S.Ls1,M.S(pt).Ls1),...
                            isequal(ptS.S.Ls2,M.S(pt).Ls2),...
                            isequal(ptS.S.Ls3,M.S(pt).Ls3),...
                            isequal(ptS.S.Ls4,M.S(pt).Ls4),...
                            isequal(ptS.S.amplitude,M.S(pt).amplitude)])
                        warning(['Local stimulation parameters for ',M.patient.list{pt},' are different to the ones stored in this Lead group analysis with the same name. Please check.']);
                    end
                end
            end
        end
        
        
        try
            setappdata(handles.leadfigure,'elstruct',elstruct);
        end
    else
        M.vilist={};
        M.fclist={};
    end
end

disp('Refreshing clinical list...');
% refresh clinical list
set(handles.clinicallist,'String',M.clinical.labels);
try set(handles.clinicallist,'Value',M.ui.clinicallist); end

if get(handles.clinicallist,'Value')>length(get(handles.clinicallist,'String'))
    set(handles.clinicallist,'Value',length(get(handles.clinicallist,'String')));
end

% store everything in Model
disp('Storing everything in model...');
if ~isempty(M.patient.list)
    t=datetime('now');
    t.Format='uuuMMddHHmmss';
    M.ui.lastupdated=str2double(char(t));
    setappdata(handles.leadfigure,'M',M);
end

disp('Done.');

ea_busyaction('off',handles.leadfigure,'group');


function S=ea_checkS(M,S,options,handles) % helper to check that S has equally many entries as M.patient.list
if ~(length(S)==length(M.patient.list))
    if isempty(S)    % will init a blank S struct
        S=ea_initializeS(['gs_',M.guid],options,handles);
        S(1:length(M.patient.list))=S(1);
    else
        %ea_error('Stimulation parameter struct not matching patient list. Lead group file potentially corrupted.');
    end
end


function ea_stats=ea_rmssstimulations(ea_stats,M)
% function that will remove all stimulations not labeled 'gs'
todel=[];
if isfield(ea_stats, 'stimulation')
    for s=1:length(ea_stats.stimulation)
        if ~strcmp(ea_stats.stimulation(s).label,['gs_',M.guid])
            todel=[todel,s];
        end
    end
    ea_stats.stimulation(todel)=[];
end
