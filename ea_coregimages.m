function affinefile = ea_coregimages(options,moving,fixed,ofile,otherfiles,writeoutmat,masks,interp)
% Generic function for image coregistration
%
% 1. Moving image will keep untouched unless the output points to the same
%    image.
% 2. 'otherfiles' will always be overwritten. If you don't want them to be
%    overwritten, use 'ea_apply_coregistration' afterwards to apply the
%    returned forward transform ('affinefile{1}') to 'otherfiles' rather
%    than supply 'otherfiles' parameter here.


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

affinefile = {};

if ~exist('masks','var') || isempty(masks)
    masks={};
end

% Only effective for SPM, 4th degree BSpline by default.
% ANTs and BRAINSFit use Linear interpolation, FSL use Sinc interpolation.
if ~exist('interp','var') || isempty(interp)
    interp = 4;
end

switch lower(options.coregmr.method)
    case lower({'ANTs (Avants 2008)', 'ANTs'})
        affinefile = ea_ants_linear(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles,masks);
    case lower({'BRAINSFit (Johnson 2007)', 'BRAINSFit'})
        affinefile = ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case lower({'FLIRT (Jenkinson 2001 & 2002)', 'FLIRT'})
        affinefile = ea_flirt(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case lower({'FLIRT BBR (Greve and Fischl 2009)', 'FLIRT BBR', 'FLIRTBBR'})
        affinefile = ea_flirtbbr(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case lower({'Hybrid SPM & ANTs', 'HybridSPMANTs'})
        % Copy moving image to out image first, since SPM will change
        % the header of the moving image.
        copyfile(moving, ofile);
        ea_spm_coreg(options,ofile,fixed,'nmi',0,otherfiles,writeoutmat)
        affinefile = ea_ants_linear(fixed,...
            ofile,...
            ofile,writeoutmat,otherfiles);
    case lower({'Hybrid SPM & BRAINSFIT', 'HybridSPMBRAINSFIT'})
        % Copy moving image to out image first, since SPM will change
        % the header of the moving image.
        copyfile(moving, ofile);
        ea_spm_coreg(options,ofile,fixed,'nmi',0,otherfiles,writeoutmat)
        affinefile = ea_brainsfit(fixed,...
            ofile,...
            ofile,writeoutmat,otherfiles);
    case lower({'Hybrid SPM & FLIRT', 'HybridSPMFLIRT'})
        % Copy moving image to out image first, since SPM will change
        % the header of the moving image.
        copyfile(moving, ofile);
        ea_spm_coreg(options,ofile,fixed,'nmi',0,otherfiles,writeoutmat)
        affinefile = ea_flirt(fixed,...
            ofile,...
            ofile,writeoutmat,otherfiles);
    case lower({'SPM (Friston 2007)', 'SPM'})
        affinefile = ea_spm_coreg(options,moving,fixed,'nmi',1,otherfiles,writeoutmat,interp);

        [fpth, fname, ext] = ea_niifileparts(moving);
        spmoutput = [fileparts(fpth), filesep, 'r', fname, ext];
        if ~strcmp(spmoutput, ofile)
            movefile(spmoutput, ofile);
        end

        for ofi=1:length(otherfiles)
            [fpth, fname, ext] = ea_niifileparts(otherfiles{ofi});
            spmoutput = [fileparts(fpth), filesep, 'r', fname, ext];
            movefile(spmoutput, [fpth, ext]);
        end
    case lower({'ANTs Nonlinear Coregistration', 'ANTsNonLinear'})
        transforms = ea_ants_nonlinear_coreg(fixed, moving, ofile, ...
            options.prefs.machine.normsettings, 'NULL', 'NULL', 'ea_antspreset_ants_wiki');
        ea_delete(transforms);
    otherwise
        warning('Coregistrion method not recognized...');
        return;
end

% Reslice images if needed (fix qform/sform issues)
V1 = ea_open_vol(fixed);
V2 = ea_open_vol(ofile);
if ~isequal(V1.mat, V2.mat)
    ea_conformspaceto(fixed, ofile, 1);
end
