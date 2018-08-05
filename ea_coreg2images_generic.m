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
        affinefile = ea_spm_coreg(options,moving,fixed,'nmi',1,otherfiles,writeoutmat,interp);
        try % will fail if ofile is same string as r mfilen..
            [pth, fn, ext] = fileparts(moving);
            movefile(fullfile(pth, ['r', fn, ext]), ofile);
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
        ea_spm_coreg(options,moving,fixed,'nmi',0,otherfiles,writeoutmat)
        affinefile = ea_ants(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Hybrid SPM & FSL' % Hybrid SPM -> FSL
        ea_spm_coreg(options,moving,fixed,'nmi',0,otherfiles,writeoutmat)
        affinefile = ea_flirt(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Hybrid SPM & BRAINSFIT' % Hybrid SPM -> Brainsfit
        ea_spm_coreg(options,moving,fixed,'nmi',0,otherfiles,writeoutmat)
        affinefile = ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
end

% ea_conformspaceto(fixed, ofile); % fix qform/sform issues.
