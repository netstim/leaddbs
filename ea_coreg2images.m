function affinefile = ea_coreg2images(options,moving,fixed,ofile,otherfiles,writeoutmat,msks,interp)

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

[directory, mfilen, ext] = ea_niifileparts(moving);
directory = [fileparts(directory),filesep];
mfilen = [mfilen,ext];

if exist([directory,'raw_',mfilen], 'file')
    copyfile([directory,'raw_',mfilen], [directory,mfilen]);
else
    copyfile([directory,mfilen], [directory,'raw_',mfilen]);
end

switch options.coregmr.method
    case 'SPM' % SPM
        affinefile = ea_spm_coreg(options,moving,fixed,'nmi',1,otherfiles,writeoutmat,interp);

        try % will fail if ofile is same string as r mfilen..
            movefile([directory,'r',mfilen],ofile);
        end

        for ofi=1:length(otherfiles)
            [pth,fn,ext] = fileparts(otherfiles{ofi});
            movefile(fullfile(pth,['r',fn,ext]), fullfile(pth,[fn,ext]));
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

ea_conformspaceto(fixed, ofile); % fix qform/sform issues.
