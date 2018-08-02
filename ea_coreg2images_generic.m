function affinefile = ea_coreg2images_generic(options,moving,fixed,ofile,otherfiles,writeoutmat,msks,interp)
% generic version of coreg2images which does not copy 'raw_' files but
% leaves original versions untouched.
if ~exist('otherfiles','var')
    otherfiles = {};
elseif isempty(otherfiles) % [] or {} or ''
    otherfiles = {};
elseif ischar(otherfiles) % single file, make it to cell string
    otherfiles = {otherfiles};
end

if ~exist('writeoutmat','var') || isempty(writeoutmat)
    writeoutmat = 0;
end

if ~exist('affinefile','var') || isempty(affinefile)
    affinefile = {};
end

if ~exist('msks','var') || isempty(msks)
    msks={};
end


if ~exist('interp','var') || isempty(interp)
    interp=4;
end

switch options.coregmr.method
    case 'SPM' % SPM
        % for SPM process everything in tmp dir.
        tmpdir=ea_getleadtempdir;
        uuid=ea_generate_uuid;
        [movingbase,movingorig]=fileparts(moving);
        [fixedbase,fixedorig]=fileparts(fixed);
        copyfile(moving,[tmpdir,uuid,'.nii']);
        moving=[tmpdir,uuid,'.nii'];

        [directory,mfilen,ext]=fileparts(moving);
        directory=[directory,filesep];
        mfilen=[mfilen,ext];


        for ofi=1:length(otherfiles)
            ofiuid{ofi}=ea_generate_uuid;
            copyfile(otherfiles{ofi},[tmpdir,ofiuid{ofi},'.nii']);
            copiedotherfiles{ofi}=[tmpdir,ofiuid{ofi},'.nii'];
        end

        if exist('copiedotherfiles','var')
            commaoneotherfiles=prepforspm(copiedotherfiles);
        else
            commaoneotherfiles={};
        end

        affinefile = ea_spm_coreg(options,appendcommaone(moving),appendcommaone(fixed),'nmi',1,commaoneotherfiles,writeoutmat,interp);
        if exist(fullfile(tmpdir,[ea_stripex(moving),'2',fixedorig,'_spm.mat']),'file')
            movefile(fullfile(tmpdir,[ea_stripex(moving),'2',fixedorig,'_spm.mat']),...
                fullfile(movingbase,[movingorig,'2',fixedorig,'_spm.mat']));
        end

        if exist(fullfile(tmpdir,[fixedorig,'2',ea_stripex(moving),'_spm.mat']),'file')
            movefile(fullfile(tmpdir,[fixedorig,'2',ea_stripex(moving),'_spm.mat']),...
                fullfile(movingbase,[fixedorig,'2',movingorig,'_spm.mat']));
        end

        try % will fail if ofile is same string as r mfilen..
            movefile([tmpdir,'r',uuid,'.nii'],ofile);
        end

        for ofi=1:length(otherfiles)
            [pth,fn,ext]=fileparts(otherfiles{ofi});
            movefile([tmpdir,'r',ofiuid{ofi},'.nii'],fullfile(pth,['r',fn,ext]));
        end

    case 'FSL' % FSL
        affinefile = ea_flirt(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'ANTs' % ANTs
        affinefile = ea_ants(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles,msks);

    case 'BRAINSFIT' % BRAINSFit
        affinefile = ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Hybrid SPM & ANTs' % Hybrid SPM -> ANTs
        commaoneotherfiles=prepforspm(otherfiles);
        ea_spm_coreg(options,appendcommaone(moving),appendcommaone(fixed),'nmi',0,commaoneotherfiles,writeoutmat)
        affinefile = ea_ants(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Hybrid SPM & FSL' % Hybrid SPM -> FSL
        commaoneotherfiles=prepforspm(otherfiles);
        ea_spm_coreg(options,appendcommaone(moving),appendcommaone(fixed),'nmi',0,commaoneotherfiles,writeoutmat)
        affinefile = ea_flirt(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Hybrid SPM & BRAINSFIT' % Hybrid SPM -> Brainsfit
        commaoneotherfiles=prepforspm(otherfiles);
        ea_spm_coreg(options,appendcommaone(moving),appendcommaone(fixed),'nmi',0,commaoneotherfiles,writeoutmat)
        affinefile = ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
end

% ea_conformspaceto(fixed, ofile); % fix qform/sform issues.


function otherfiles=prepforspm(otherfiles)
if size(otherfiles,1)<size(otherfiles,2)
    otherfiles=otherfiles';
end

for fi=1:length(otherfiles)
    otherfiles{fi}=appendcommaone(otherfiles{fi});
end


function fname=appendcommaone(fname)
if ~strcmp(fname(end-1:end),',1')
    fname= [fname,',1'];
end
