function affinefile = ea_coreg2images(options,moving,fixed,ofile,otherfiles,writeoutmat,msks)

if nargin < 5
    otherfiles = {};
elseif isempty(otherfiles) % [] or {} or ''
    otherfiles = {};
elseif ischar(otherfiles) % single file, make it to cell string
    otherfiles = {otherfiles};
end

if nargin < 6
    writeoutmat = 0;
    affinefile = {};
end
if ~exist('msks','var')
    msks={};
end
[directory,mfilen,ext]=fileparts(moving);
directory=[directory,filesep];
mfilen=[mfilen,ext];

if exist([directory,'raw_',mfilen],'file')
    copyfile([directory,'raw_',mfilen],[directory,mfilen]);
else
    copyfile([directory,mfilen],[directory,'raw_',mfilen]);
end

switch options.coregmr.method
    case 'SPM' % SPM
        commaoneotherfiles=prepforspm(otherfiles);

        affinefile = ea_docoreg_spm(options,appendcommaone(moving),appendcommaone(fixed),'nmi',1,commaoneotherfiles,writeoutmat);
        try % will fail if ofile is same string as r mfilen..
            movefile([directory,'r',mfilen],ofile);
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
