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

tmpdir=ea_getleadtempdir;
uid=ea_generate_guid;
copyfile(moving,[tmpdir,uid,'.nii']);
moving=[tmpdir,uid,'.nii'];

[directory,mfilen,ext]=fileparts(moving);
directory=[directory,filesep];
mfilen=[mfilen,ext];


if ~exist('interp','var') || isempty(interp)
    interp=4;
end

switch options.coregmr.method
    case 'SPM' % SPM
        commaoneotherfiles=prepforspm(otherfiles);

        affinefile = ea_docoreg_spm(options,appendcommaone(moving),appendcommaone(fixed),'nmi',1,commaoneotherfiles,writeoutmat,interp);
        try % will fail if ofile is same string as r mfilen..
            movefile([tmpdir,'r',uid,'.nii'],ofile);
        end
        for ofi=1:length(otherfiles)
            [pth,fn,ext]=fileparts(otherfiles{ofi});
            movefile(fullfile(pth,['r',fn,ext]),fullfile(pth,[fn,ext]));
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
        ea_docoreg_spm(options,appendcommaone(moving),appendcommaone(fixed),'nmi',0,commaoneotherfiles,writeoutmat)
        affinefile = ea_ants(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Hybrid SPM & FSL' % Hybrid SPM -> FSL
        commaoneotherfiles=prepforspm(otherfiles);
        ea_docoreg_spm(options,appendcommaone(moving),appendcommaone(fixed),'nmi',0,commaoneotherfiles,writeoutmat)
        affinefile = ea_flirt(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Hybrid SPM & BRAINSFIT' % Hybrid SPM -> Brainsfit
        commaoneotherfiles=prepforspm(otherfiles);
        ea_docoreg_spm(options,appendcommaone(moving),appendcommaone(fixed),'nmi',0,commaoneotherfiles,writeoutmat)
        affinefile = ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
end
ea_conformspaceto(fixed, ofile); % fix qform/sform issues.


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
