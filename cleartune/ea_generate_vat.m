function [Efields,S]=ea_generate_optim_vat(varargin)

patselect = varargin{1};
ampselect = varargin{2};
constcurr = varargin{4};
side = varargin{5};
writeVTA = varargin{6};
modelVTA = varargin{7};
%va = constcurr;
%concval = varargin{3}(1:8);
concval = varargin{3}(1:8);
%concval = [concval,100]; %switch this off for bipolar
%% Load and define options
options = ea_setopts_local;
options.native = 0;
options.groupmode = 1;
options.groupid = 'cleartune';

va = constcurr; % 0 for constant curr
                % 1 for voltage

% amps.min = min_amplitude;
% amps.max = max_amplitude;
% amps.stepsize = stepsize;
% multset = 'mon';

resultfig=figure('visible','off');

%% Start iterating through patients

options = ea_getptopts(patselect, options);

fprintf('\nProcessing sub-%s...\n\n', options.subj.subjId);

[coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(options);
elstruct(1).coords_mm=coords_mm;
elstruct(1).coords_acpc=coords_acpc;
elstruct(1).trajectory=trajectory;
elstruct(1).name = ['sub-', options.subj.subjId];
elstruct(1).markers=markers;

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
options.patientname = options.subj.subjId;
setappdata(resultfig,'elstruct',elstruct);
setappdata(resultfig,'options',options);
setappdata(resultfig,'elspec',options.elspec);
setappdata(resultfig,'resultfig',resultfig);

%% Define stimulation settings

%allstims{1} = repmat(((1:ncnttmp)-1)',1,size(amps.min:amps.stepsize:amps.max,2));
%allstims{2} = repmat(amps.min:amps.stepsize:amps.max,ncnttmp,1);

%runs = 1:numel(allstims{1});
t=load([ea_getearoot,'templates',filesep,'electrode_models',filesep,options.elspec.matfname '.mat']); % defines electrode variable
elt=load([ea_getearoot,'templates',filesep,'standard_efields' filesep 'standard_efield_' options.elspec.matfname '.mat']);
whichContact = find(concval);
if length(whichContact) > 1
    whichContact = [num2str(whichContact(1:end-1))];
    whichContact = whichContact(~isspace(whichContact));
end
if side == 1
    stimtmpR = [ampselect,concval];
    tmp_ampsel = 0;
    tmp_concval = zeros(1,length(concval));
    stimtmpL = [tmp_ampsel,tmp_concval];
elseif side == 2
    stimtmpL = [ampselect,concval];
    tmp_ampsel = 0;
    tmp_concval = zeros(1,length(concval));
    stimtmpR = [tmp_ampsel,tmp_concval];
end
    
S = ea_initializeS(options);
S = ea_cleartune_generateMfile(stimtmpR,stimtmpL,S,va);
S.label = ['amp_R_L_',num2str(ampselect,'%.2f'),'_',num2str(ampselect,'%.2f'),'_contactR_L_',num2str(whichContact,'%d'),'_',num2str(whichContact,'%d')];

% Define the name of the folder for the nii to be saved in
if va == 0
    volcur = 'mA';
elseif va == 1
    volcur = 'V';
else
    volcur = '??';
end
fname = [volcur, '_', num2str(round(options.prefs.machine.vatsettings.horn_cgm*100),'%02d'), '_', num2str(round(options.prefs.machine.vatsettings.horn_cwm*100),'%02d')];

for hem=side
    if hem == side
        disp([' ', newline, 'Patient ', patselect, newline,'Simulating efield: ', fname, ' side ', num2str(hem),' | ', S.label])
    end
    setappdata(resultfig,'elstruct',elstruct(1));
    setappdata(resultfig,'elspec',options.elspec);
    if strcmp(modelVTA,'Fastfield')
        Efields=ea_genvat_cleartune_fastfield(S,hem,options,fname,resultfig,t.electrode,elt);
        
        if hem == side && writeVTA
            ea_write_nii(Efields)
            
        end
    else
        tic;
        outputEfield = ea_genvat_cleartune_horn('',S,hem,options,fname,resultfig);
        if isempty(outputEfield)
            Vvate = createEmptyNii(t,elstruct,elt,side,fname);
            Efields= Vvate;
        else
            Efields = outputEfield;
        end
        toc;
        return % one E-field at a time
    end
end


close(resultfig);
function Vvate = createEmptyNii(t,elstruct,elt,side,fname)
    % create nifti
    [~, ~, endian] = computer;
    switch endian
        case 'L'
            endian = 0;
        case 'B'
            endian = 1;
    end
    [trans_mat] = get_trans_mat(t.electrode,elstruct,elt.grid_vec,side);
    gv=elt.grid_vec;
    Vvate.img = zeros(100,100,100);
    res=100;
    chun1=randperm(res); chun2=randperm(res); chun3=randperm(res);
    Vvate.mat=mldivide([(chun1);(chun2);(chun3);ones(1,res)]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,res)]')';
    Vvate.mat = trans_mat * Vvate.mat;
    Vvate.dim=[res,res,res];
    Vvate.dt = [16, endian];
    Vvate.n=[1 1];
    Vvate.fname = fname;
    Vvate.descrip='lead dbs - vat';
return


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