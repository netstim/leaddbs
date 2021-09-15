function varargout=ea_assignbackdrop(bdstring,options,subpat,native)

if ~exist('subpat','var')
    subpat='Patient';
end

if ~exist('native','var')
    if isfield(options,'native')
        native = options.native;
    else
        native = 0;   % default
    end
end

BDlist=getbdlist;

if strcmp(bdstring, 'list')
    % determine whether we are in No patient mode (could be called from
    % lead group or called from an empty patient viewer / lead anatomy
    if ~exist('options','var')
        options.patientname='';
    end

    % by default assuming patient mode, preop and postop images exist.
    nopatientmode=0;
    haspreop=1;
    haspostop=1;

    % check patient mode
    if isfield(options,'groupmode')
        nopatientmode=options.groupmode;
    elseif ~isfield(options,'patientname')
        nopatientmode=1;
    elseif strcmp(options.patientname,'No Patient Selected')
        nopatientmode=1;
    end

    % check if preop and postop images exist
    if ~nopatientmode
        if isempty(dir([options.root,options.patientname,filesep,'anat_*.nii']))
            haspreop=0;
        end
        try
            % use this as a probe to see if required patient postop images exist.
            assignpatspecific(options, native);
        catch
            haspostop=0;
        end
    end

    if ~haspostop && ~haspreop
        nopatientmode=1;
    end

    if nopatientmode
        varargout{1}=[ea_standardspacelist];
    else
        if haspreop
            [~, preopfiles]=ea_assignpretra(options);
            preop=cellfun(@(x) [subpat, ' Pre-OP (',upper(regexp(x, '(?<=anat_)(.*)(?=\.nii)', 'match', 'once')), ')'], preopfiles, 'Uniform', 0)';
        else
            preop={''};
        end
        postop = {[subpat,' Post-OP']};
        if native
            varargout{1}=[preop,...
                postop(logical(haspostop))];
        else
            varargout{1}=[ea_standardspacelist,...
                preop,...
                postop(logical(haspostop))];
        end

    end

    if ~native
        % check for additional template backdrops
            for bd=1:length(BDlist{1})
                varargout{1}=[varargout{1},...
                    BDlist{2}{bd}];
            end
    end

    % add manual choose:
    varargout{1}=[varargout{1},...
                {'Choose...'}];

elseif regexp(bdstring, ['^', subpat,' Pre-OP \(.*\)$'])    % pattern: "Patient Pre-OP (*)"
    whichpreop=lower(regexp(bdstring, ['(?<=^', subpat,' Pre-OP \()(.*)(?=\))'],'match','once'));
    options=ea_tempswitchoptstopre(options, native, whichpreop);
    [Vtra,Vcor,Vsag]=assignpatspecific(options, native);
    varargout{1}=Vtra;
    varargout{2}=Vcor;
    varargout{3}=Vsag;

elseif strcmp(bdstring, [subpat, ' Post-OP'])
    [Vtra,Vcor,Vsag]=assignpatspecific(options, native);
    varargout{1}=Vtra;
    varargout{2}=Vcor;
    varargout{3}=Vsag;

elseif strcmp(bdstring, 'BigBrain 100 um ICBM 152 2009b Sym (Amunts 2013)')
%     if ~ea_checkinstall('bigbrain',0,1)
%         ea_error('BigBrain is not installed and could not be installed automatically. Please make sure that Matlab is connected to the internet.');
%     end
%     varargout{1}=spm_vol(fullfile(ea_space(options),'bigbrain_2015_100um_bb.nii'));
%     varargout{2}=varargout{1};
%     varargout{3}=varargout{1};

elseif regexp(bdstring, ['^',ea_getspace,' '])    % pattern: "MNI_ICBM_2009b_NLIN_ASYM *"
    spacedef=ea_getspacedef;
    template=lower(strrep(strrep(bdstring,[ea_getspace,' '],''),[' (',spacedef.citation{1},')'],''));

    varargout{1}=spm_vol(ea_niigz(fullfile(ea_space(options),[template])));
    varargout{2}=varargout{1};
    varargout{3}=varargout{1};
elseif strcmp(bdstring,'Choose...')
    keyboard

elseif ismember(bdstring,BDlist{2})
    [~,ix]=ismember(bdstring,BDlist{2});
    varargout{1}=ea_load_nii([ea_space,'backdrops',filesep,BDlist{1}{ix}]);
    varargout{2}=varargout{1};
    varargout{3}=varargout{1};
else    % custom backdrop file
    varargout{1}=spm_vol(bdstring);
    varargout{2}=varargout{1};
    varargout{3}=varargout{1};
end


function [Vtra,Vcor,Vsag]=assignpatspecific(options, native)
if native
    switch options.modality
        case 1 % MR
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii_unnormalized));
            if exist(fullfile(options.root,options.prefs.patientdir,options.prefs.cornii_unnormalized), 'file')
                Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.cornii_unnormalized));
            else
                Vcor=Vtra;
            end
            if exist(fullfile(options.root,options.prefs.patientdir,options.prefs.sagnii_unnormalized),'file')
                Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.sagnii_unnormalized));
            else
                Vsag=Vtra;
            end
        case 2 % CT
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tp_ctnii_coregistered));
            Vcor=Vtra;
            Vsag=Vtra;
    end
else
    switch options.modality
        case 1 % MR
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gtranii));
            if exist(fullfile(options.root,options.prefs.patientdir,options.prefs.gcornii),'file')
                Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gcornii));
            else
                Vcor=Vtra;
            end
            if exist(fullfile(options.root,options.prefs.patientdir,options.prefs.gsagnii),'file')
                Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gsagnii));
            else
                Vsag=Vtra;
            end
        case 2 % CT
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tp_gctnii));
            Vcor=Vtra;
            Vsag=Vtra;
    end
end


function standardlist=ea_standardspacelist
spacedef=ea_getspacedef;
standardlist=cell(1,length(spacedef.templates));
for t=1:length(spacedef.templates)
	standardlist{t}=[spacedef.name,' ',upper(spacedef.templates{t}),' (',spacedef.citation{1},')'];
end


function BDlist=getbdlist
BDlist=[{[]},{[]}]; % empty.
if exist([ea_space,'backdrops',filesep,'backdrops.txt'],'file')
    fid=fopen([ea_space,'backdrops',filesep,'backdrops.txt']);
    BDlist=textscan(fid,'%s %s');
    BDlist{2}=ea_underscore2space(BDlist{2});
end


function options=ea_tempswitchoptstopre(options, native, whichpreop)
% this generates a very temporary fake options struct that points to preop
% data instead of postop data.

if native
    options.prefs.tranii_unnormalized=['anat_',whichpreop,'.nii'];
    options.prefs.cornii_unnormalized=['anat_',whichpreop,'.nii'];
    options.prefs.sagnii_unnormalized=['anat_',whichpreop,'.nii'];
    options.prefs.tp_ctnii_coregistered=['anat_',whichpreop,'.nii'];
else
    [options,preniis]=ea_assignpretra(options);
    gfi=['glanat','.nii']; % default use the file available.
    try
        if strcmp(['anat_',whichpreop,'.nii'], preniis{1})  % whichpreop: "t1", preniis{1}: "anat_t1.nii"
            gfi=['glanat','.nii'];
        else
            gfi=['glanat_',whichpreop,'.nii'];
        end
    end
    options.prefs.gtranii=gfi;
    options.prefs.gcornii=gfi;
    options.prefs.gsagnii=gfi;
    options.prefs.tp_gctnii=gfi;

    if ~exist([options.root,options.patientname,filesep,gfi],'file')
        to{1}=[options.root,options.patientname,filesep,gfi];
        from{1}=[options.root,options.patientname,filesep,['anat_',whichpreop,'.nii']];
        ea_apply_normalization_tofile(options,from,to,0);
    end
end
