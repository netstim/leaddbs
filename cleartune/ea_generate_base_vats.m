function [Efields,S]=ea_generate_base_vats(varargin)

patselect=varargin{1};
min_amplitude=varargin{2};
max_amplitude=varargin{3};
stepsize=varargin{4};
va=varargin{5};

%% Load and define options
options = ea_setopts_local;
options.native = 0;
options.groupmode = 1;
options.groupid = 'cleartune';

switch va
    case 'mA'
        va = 0;             % Constant current stimulation
    case 'V'
        va = 1;             % Constant voltage stimulation
end

amps.min = min_amplitude;
amps.max = max_amplitude;
amps.stepsize = stepsize;
multset = 'mon';

resultfig=figure('visible','off');

%% Start iterating through patients
for pt=1:length(patselect)
    options = ea_getptopts(patselect{pt}, options);

    fprintf('\nProcessing sub-%s...\n\n', options.subj.subjId);

    [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(options);
    elstruct(pt).coords_mm=coords_mm;
    elstruct(pt).coords_acpc=coords_acpc;
    elstruct(pt).trajectory=trajectory;
    elstruct(pt).name = ['sub-', options.subj.subjId];
    elstruct(pt).markers=markers;

    options.numcontacts=size(coords_mm{1},1);
    options.d3.verbose='off';
    options.d3.elrendering=1;	% hard code to viz electrodes in this setting.
    options.d3.exportBB=0;	% don't export brainbrowser struct by default
    options.d3.colorpointcloud=0;
    options.d3.hlactivecontacts=1;
    options.d3.showactivecontacts =1;
    options.d3.showpassivecontacts=1;
    options.d3.exportBB=0;
    options.expstatvat.do=0;
    options.leadprod = 'group';
    options.patient_list=patselect;
    options.d3.mirrorsides=0;
    options.atlasset = options.prefs.machine.vatsettings.horn_atlasset;

    setappdata(resultfig,'elstruct',elstruct);
    setappdata(resultfig,'options',options);
    setappdata(resultfig,'elspec',options.elspec);
    setappdata(resultfig,'resultfig',resultfig);

    %% Define stimulation settings
    ncnttmp = options.numcontacts;

    allstims{1} = repmat(((1:ncnttmp)-1)',1,size(amps.min:amps.stepsize:amps.max,2));
    allstims{2} = repmat(amps.min:amps.stepsize:amps.max,ncnttmp,1);

    runs = 1:numel(allstims{1});
    t=load([ea_getearoot,'templates',filesep,'electrode_models',filesep,options.elspec.matfname '.mat']); % defines electrode variable
    elt=load([ea_getearoot,'templates',filesep,'standard_efields' filesep 'standard_efield_' options.elspec.matfname '.mat']);

    %% Iterate through all stimulations
    for run = runs
            stimtmp = zeros(1,9);
            stimtmp(1,1) = allstims{1,2}(run);
            stimtmp(1,allstims{1,1}(run)+2) = -100;
            options.UsePreExistingStim = 0;
            S(run) = ea_initializeS(options);
            S(run) = ea_cleartune_generateMfile(stimtmp,stimtmp,S(run),va);
            S(run).label = ['c',num2str(allstims{1}(run),'%02d'),'_a',num2str(allstims{2}(run)*10,'%02d')];

        % Define the name of the folder for the nii to be saved in
        if va == 0
            volcur = 'mA';
        elseif va == 1
            volcur = 'V';
        else
            volcur = '??';
        end

        fname = [volcur, '_', num2str(round(options.prefs.machine.vatsettings.horn_cgm*100),'%02d'), '_', num2str(round(options.prefs.machine.vatsettings.horn_cwm*100),'%02d')];
        %% Iterate both hemispheres
        for side=1:2
            disp([' ', newline, 'Patient ', patselect{pt}, newline,'Simulating efield: ', fname, ' side ', num2str(side),' | ', S(run).label])
            setappdata(resultfig,'elstruct',elstruct(pt));
            setappdata(resultfig,'elspec',options.elspec);
            Efields(run,side)=ea_genvat_cleartune_fastfield(S(run),side,options,fname,resultfig,t.electrode,elt);
        end
    end
end

close(resultfig);


function options=ea_setopts_local

options.earoot=ea_getearoot;
options.verbose=3;
options.sides=1:2; % re-check this later..
options.fiberthresh=1;
options.writeoutstats=1;
options.writeoutpm = 0;
options.colormap=jet;
options.d3.write=1;
options.d3.prolong_electrode=2;
options.d3.writeatlases=1;
options.macaquemodus=0;

% try
%    options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
% catch % too many entries..
%    set(handles.atlassetpopup,'Value',1);
%    options.atlasset=1;
% end

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
