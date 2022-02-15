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
        nopatientmode = options.groupmode;
    elseif isempty(options.subj.subjId)
        nopatientmode=1;
    end

    % check if preop and postop images exist
    if ~nopatientmode
        if ~isfile(options.subj.coreg.anat.preop.(options.subj.AnchorModality))
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
        varargout{1}=ea_standardspacelist;
    else
        if haspreop
            preopfiles = fieldnames(options.subj.coreg.anat.preop);
            preop=cellfun(@(x) [subpat, ' Pre-OP (' x ')'], preopfiles, 'Uniform', 0)';
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
    whichpreop = (regexp(bdstring, ['(?<=^', subpat,' Pre-OP \()(.*)(?=\))'],'match','once'));
    options = ea_switchpost2pre(options, native, whichpreop);
    [Vtra,Vcor,Vsag] = assignpatspecific(options, native);
    varargout{1} = Vtra;
    varargout{2} = Vcor;
    varargout{3} = Vsag;

elseif strcmp(bdstring, [subpat, ' Post-OP'])
    [Vtra,Vcor,Vsag] = assignpatspecific(options, native);
    varargout{1} = Vtra;
    varargout{2} = Vcor;
    varargout{3} = Vsag;

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

    varargout{1}=spm_vol(ea_niigz(fullfile(ea_space(options),template)));
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


function [Vtra,Vcor,Vsag] = assignpatspecific(options, native)
if native
    switch options.modality
        case 1 % MR
            Vtra = spm_vol(options.subj.coreg.anat.postop.ax_MRI);
            if isfield(options.subj.coreg.anat.postop, 'cor_MRI') && ifile(options.subj.coreg.anat.postop.cor_MRI)
                Vcor = spm_vol(options.subj.coreg.anat.postop.cor_MRI);
            else
                Vcor = Vtra;
            end
            if isfield(options.subj.coreg.anat.postop, 'sag_MRI') && ifile(options.subj.coreg.anat.postop.sag_MRI)
                Vsag = spm_vol(options.subj.coreg.anat.postop.sag_MRI);
            else
                Vsag = Vtra;
            end
        case 2 % CT
            Vtra = spm_vol(options.subj.coreg.anat.postop.tonemapCT);
            Vcor = Vtra;
            Vsag = Vtra;
    end
else
    switch options.modality
        case 1 % MR
            Vtra = spm_vol(options.subj.norm.anat.postop.ax_MRI);
            if isfield(options.subj.norm.anat.postop, 'cor_MRI') && ifile(options.subj.norm.anat.postop.cor_MRI)
                Vcor = spm_vol(options.subj.norm.anat.postop.cor_MRI);
            else
                Vcor = Vtra;
            end
            if isfield(options.subj.norm.anat.postop, 'sag_MRI') && ifile(options.subj.norm.anat.postop.sag_MRI)
                Vsag = spm_vol(options.subj.norm.anat.postop.sag_MRI);
            else
                Vsag = Vtra;
            end
        case 2 % CT
            Vtra = spm_vol(options.subj.norm.anat.postop.tonemapCT);
            Vcor = Vtra;
            Vsag = Vtra;
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


function options = ea_switchpost2pre(options, native, whichpreop)
% Generates a very temporary fake options struct that links pre-op to
% post-op data.

if native
    options.subj.coreg.anat.postop.ax_MRI = options.subj.coreg.anat.preop.(whichpreop);
    options.subj.coreg.anat.postop.cor_MRI = options.subj.coreg.anat.preop.(whichpreop);
    options.subj.coreg.anat.postop.sag_MRI = options.subj.coreg.anat.preop.(whichpreop);
    options.subj.coreg.anat.postop.tonemapCT = options.subj.coreg.anat.preop.(whichpreop);
else
    normImage = options.subj.preopAnat.(options.subj.AnchorModality).norm;
    normImage = strrep(normImage, ['_', options.subj.AnchorModality, '.'], ['_', whichpreop, '.']);

    if ~isfile(normImage)
        to{1} = normImage;
        from{1} = options.subj.preopAnat.(whichpreop).coreg;
        ea_apply_normalization_tofile(options,from,to,0);
    end

    options.subj.norm.anat.postop.ax_MRI = normImage;
    options.subj.norm.anat.postop.cor_MRI = normImage;
    options.subj.norm.anat.postop.sag_MRI = normImage;
    options.subj.norm.anat.postop.tonemapCT = normImage;
end
