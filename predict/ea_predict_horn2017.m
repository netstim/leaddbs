function varargout=ea_predict_horn2017(varargin)

specs.modelname='Horn et al., 2017 Annals of Neurology';
specs.modelshortname='horn2017';
specs.default.dMRIcon='HCP_MGH_30fold_groupconnectome (Horn 2017)';
specs.default.fMRIcon='GSP 1000 (Yeo 2011)>Full Set (Yeo 2011)';
specs.feats={'dMRI','fMRI'}; % could be Coords and VTA
specs.metrics={'% UPDRS-III Improvement'};
specs.support={'HCP_MGH_30fold_groupconnectome (Horn 2017)','PPMI_90 (Ewert 2017)','PPMI 74_15 (Horn 2017)'};

if nargin==1 && ischar(varargin{1}) && strcmp(varargin{1},'specs')
    varargout{1}=specs;
    return
end

options=varargin{1};
pt=varargin{2}; % patient number (of options.uivatdirs) defined in outer loop.
directory=[options.uivatdirs{pt},filesep];
load(fullfile(ea_getearoot,'predict','models','horn2017_AoN','modeldata.mat'));

[~, subPrefix] = fileparts([options.uivatdirs{pt}, '_']);
fConnName = ea_getConnLabel(options.predict.fMRIcon);
dConnName = ea_getConnLabel(options.predict.dMRIcon);
fMRIMapName = [subPrefix, 'sim-binary_model-simbio_seed-fMRI_conn-', fConnName, '_desc-AvgRFz_funcmap.nii'];
dMRIMapName = [subPrefix, 'sim-binary_model-simbio_seed-dMRI_conn-', dConnName, '_strucmap.nii'];
SKdMRIMapName = [subPrefix, 'sim-binary_model-simbio_seed-dMRI_conn-', dConnName, '_desc-NormSmooth_strucmap.nii']; % Smoothed and normalized

feats=[0,0];
stimname=options.predict.stimulation;

%% get seed maps of VTAs
if ismember('dMRI',options.predict.includes)
    feats(1)=1;
    if startsWith(options.predict.dMRIcon,'Precomputed: ')
        options.predict.dMRIcon = erase(options.predict.dMRIcon, 'Precomputed: ');
    else
        % -> run connectome mapper on patient
        run_mapper_vat_local(uivatdirs{pt},stimname,0,options.predict.dMRIcon,1,options.predict.fMRIcon)
    end
    if ~exist([options.uivatdirs{pt},filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,SKdMRIMapName],'file')
        ea_dosk([options.uivatdirs{pt},filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,dMRIMapName],modeldata.mask)
    end
    dMRImap=ea_load_nii([options.uivatdirs{pt},filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,SKdMRIMapName]);
end

if ismember('fMRI',options.predict.includes)
    feats(2)=1;
    if startsWith(options.predict.fMRIcon,'Precomputed: ')
        options.predict.fMRIcon = erase(options.predict.fMRIcon, 'Precomputed: ');
    else
        % -> run connectome mapper on patient
        run_mapper_vat_local(uivatdirs{pt},stimname,1,options.predict.dMRIcon,0,options.predict.fMRIcon)
    end
    fMRImap=ea_load_nii([options.uivatdirs{pt},filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,fMRIMapName]);
end

% model
warn=0;
if ~exist(fullfile(ea_getearoot,'predict','models','horn2017_AoN','combined_maps','dMRI',options.predict.dMRIcon),'dir') && feats(1)
    dMRIcon=specs.default.dMRIcon;
    warn=1;
else
    dMRIcon=options.predict.dMRIcon;
end

if ~exist(fullfile(ea_getearoot,'predict','models','horn2017_AoN','combined_maps','fMRI',options.predict.fMRIcon),'dir') && feats(2)
    fMRIcon=specs.default.fMRIcon;
    warn=1;
else
    fMRIcon=options.predict.fMRIcon;
end

%% load canonical models and compare
if feats(1)
    dMRImod=ea_load_nii(fullfile(ea_getearoot,'predict','models','horn2017_AoN','combined_maps','dMRI',dMRIcon,'dMRI_optimal.nii'));
    dMRIsim=corr(dMRImod.img(modeldata.mask),dMRImap.img(modeldata.mask),'type','spearman','rows','pairwise');
end

if feats(2)
    fMRImod=ea_load_nii(fullfile(ea_getearoot,'predict','models','horn2017_AoN','combined_maps','fMRI',strrep(fMRIcon,'>','_'),'fMRI_optimal.nii'));
    modelvals=fMRImod.img(modeldata.mask);
    ptvals=fMRImap.img(modeldata.mask);
    infs=isinf(modelvals);
    infs=logical(infs+isinf(ptvals));
    modelvals(infs)=[];
    ptvals(infs)=[];
    fMRIsim=corr(modelvals,ptvals,'type','pearson','rows','pairwise');
end

% solve regression model
X=[modeldata.connectomes.(rmbracketspace(specs.default.dMRIcon)).dMRIsims',modeldata.connectomes.(rmbracketspace(specs.default.fMRIcon)).fMRIsims'];
if isfield(modeldata.connectomes,rmbracketspace(dMRIcon))
    X(:,1)=modeldata.connectomes.(rmbracketspace(dMRIcon)).dMRIsims';
else
    if feats(1)
        warn=1;
    end
end

if isfield(modeldata.connectomes,rmbracketspace(fMRIcon))
    X(:,2)=modeldata.connectomes.(rmbracketspace(fMRIcon)).fMRIsims';
else
    if feats(2)
        warn=1;
    end
end

if warn
    ea_warning('Please note that this prediction module is not validated for use with the selected connectome (or patient''s data). You may proceed but the results may not be meaningful.');
end

X=X(:,logical(feats));
[beta,dev,stats]=glmfit(X,modeldata.updrs3percimprov);

UpdrsHat=ea_addone(X)*beta;

avgerror=mean(abs(UpdrsHat-modeldata.updrs3percimprov));
% %R=corr(UpdrsHat,modeldata.updrs3percimprov)

Xpt=[0,0];
if feats(1)
    Xpt(1)=dMRIsim;
end

if feats(2)
    Xpt(2)=fMRIsim;
end

Xpt=Xpt(logical(feats));

updrshat=ea_addone(Xpt)*beta; % percent UPDRS-III improvement prediction

%% build improvement report:

if ~exist([directory,'predictions'],'dir')
    mkdir([directory,'predictions']);
end

mkdir([directory,'predictions',filesep,specs.modelshortname]);
hasdMRI=0; hasfMRI=0;

if ismember('dMRI',options.predict.includes)
    hasdMRI=1;
end

if ismember('fMRI',options.predict.includes)
    hasfMRI=1;
end

if hasfMRI
    cfg.fMRI.model=ea_getsurficeplots(fullfile(ea_getearoot,'predict','models','horn2017_AoN','combined_maps','fMRI',strrep(fMRIcon,'>','_'),'fMRI_optimal.nii'),[0.01,0.1,-0.01,-0.1]);
    cfg.fMRI.vta=ea_getsurficeplots(fMRImap.fname,[0.01,0.1,-0.01,-0.1]);
end

if hasdMRI
    cfg.dMRI.model=ea_getsurficeplots(fullfile(ea_getearoot,'predict','models','horn2017_AoN','combined_maps','dMRI',dMRIcon,'dMRI_optimal.nii'),[-1,0.6,nan,nan]);
    cfg.dMRI.vta=ea_getsurficeplots(dMRImap.fname,[-1,0.6,nan,nan]);
end

cfg.res.updrs3imp=updrshat;
cfg.res.updrs3err=avgerror;
load([directory,'stimulations',filesep,ea_nt(options),stimname,filesep,'stimparameters.mat'])
% cfg.stim=ea_activecontacts(S);
cfg.stim=S;
[~,cfg.stim.patientname]=fileparts(fileparts(directory));
ea_presults_horn2017(cfg);
%ea_screenshot([directory,'predictions',filesep,specs.modelshortname,filesep,'prediction.png']);


function str=rmbracketspace(str)
str=strrep(str,' ','_');
str=strrep(str,'(','_');
str=strrep(str,')','_');
str=strrep(str,'>','_');


function run_mapper_vat_local(ptdir,stimname,struc,strucc,func,funcc)
% - Lead-DBS Job created on 21-Oct-2017 19:02:11 -
% --------------------------------------
options = getoptslocal;
options.lcm.struc.do = struc;
options.lcm.func.do = func;

options.lcm.seeds = stimname;
options.lcm.struc.connectome = strucc;
options.lcm.func.connectome = funcc;
options.uivatdirs = {ptdir};
options.prefs=ea_prefs;

options.leadprod = 'mapper';

ea_run('run', options);


function options = getoptslocal
options.endtolerance = 10;
options.sprungwert = 4;
options.refinesteps = 0;
options.tra_stdfactor = 0.9;
options.cor_stdfactor = 1;
options.earoot = ea_getearoot;
options.normalize.do = 0;
options.normalize.method = [];
options.checkreg = false;
options.coregmr.check = 0;
options.coregmr.do = 0;
options.coregmr.method = '';
options.coregct.do = 0;
options.modality = 1;
options.verbose = 3;
options.sides = [1 2];
options.doreconstruction = 0;
options.autoimprove = 0;
options.axiscontrast = 8;
options.zresolution = 10;
options.atl.genpt = 0;
options.atl.can = 1;
options.atl.pt = 0;
options.atl.ptnative = 0;
options.native = 0;
options.d2.col_overlay = 1;
options.d2.con_overlay = 1;
options.d2.con_color = [1 1 1];
options.d2.lab_overlay = 1;
options.d2.bbsize = 10;
options.d2.backdrop = 'MNI152NLin2009bAsym T1 (Fonov 2011)';
options.d2.fid_overlay = 0;
options.d2.write = 0;
options.d2.atlasopacity = 0.15;
options.refinelocalization = 0;
options.scrf = 0;
options.d3.write = 0;
options.d3.prolong_electrode = 2;
options.d3.verbose = 'on';
options.d3.elrendering = 1;
options.d3.exportBB = 0;
options.d3.hlactivecontacts = 0;
options.d3.showactivecontacts = 1;
options.d3.showpassivecontacts = 1;
options.d3.showisovolume = 0;
options.d3.isovscloud = 0;
options.d3.mirrorsides = 0;
options.d3.autoserver = 0;
options.d3.expdf = 0;
optiosns.numcontacts = 4;
options.writeoutpm = 0;
options.expstatvat.do = 0;
options.fiberthresh = 10;
options.writeoutstats = 1;
options.dolc = 0;
options.lcm.seeddef = 'vats';
options.lcm.odir = '';
options.lcm.omask = [];
options.lcm.struc.espace = 1;
options.lcm.cmd = 1;
options.uipatdirs = {''};
options.lc.general.parcellation = 'AICHA reordered (Joliot 2015)';
options.lc.graph.struc_func_sim = 0;
options.lc.graph.nodal_efficiency = 0;
options.lc.graph.eigenvector_centrality = 0;
options.lc.graph.degree_centrality = 0;
options.lc.graph.fthresh = NaN;
options.lc.graph.sthresh = NaN;
options.lc.func.compute_CM = 0;
options.lc.func.compute_GM = 0;
options.lc.func.prefs.TR = 2.69;
options.lc.struc.compute_CM = 0;
options.lc.struc.compute_GM = 0;
options.lc.struc.ft.method = 'ea_ft_mesotracking_reisert';
options.lc.struc.ft.do = 1;
options.lc.struc.ft.normalize = 0;
options.exportedJob = 1;
