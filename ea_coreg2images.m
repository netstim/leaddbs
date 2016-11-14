function affinefile = ea_coreg2images(options,moving,fixed,ofile,otherfiles,writeoutmat)

if nargin < 5
    otherfiles={''}; 
elseif ischar(otherfiles)
    otherfiles = {otherfiles};
end

if nargin < 6
    writeoutmat = 0;
    affinefile = {''};
end

[directory,mfilen,ext]=fileparts(moving);
directory=[directory,filesep];
mfilen=[mfilen,ext];

if exist([directory,'raw_',mfilen],'file');
    copyfile([directory,'raw_',mfilen],[directory,mfilen]);
else
    copyfile([directory,mfilen],[directory,'raw_',mfilen]);
end

switch options.coregmr.method
    case 'Coreg MRIs: SPM' % SPM
        affinefile = ea_docoreg_spm(moving,fixed,'nmi',1,otherfiles,writeoutmat);
        movefile([directory,'r',mfilen],ofile);
        for ofi=1:length(otherfiles)
            [pth,fn,ext]=fileparts(otherfiles{ofi});
            movefile(fullfile(pth,['r',fn,ext]),fullfile(pth,[fn,ext]));
        end
    case 'Coreg MRIs: FSL' % FSL 
        affinefile = ea_flirt(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);        
    case 'Coreg MRIs: ANTs' % ANTs
        affinefile = ea_ants(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Coreg MRIs: BRAINSFIT' % BRAINSFit
        affinefile = ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Coreg MRIs: Hybrid SPM & ANTs' % Hybrid SPM -> ANTs
        ea_docoreg_spm(moving,fixed,'nmi',0,otherfiles)
        affinefile = ea_ants(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Coreg MRIs: Hybrid SPM & FSL' % Hybrid SPM -> FSL
        ea_docoreg_spm(moving,fixed,'nmi',0,otherfiles)
        affinefile = ea_flirt(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Coreg MRIs: Hybrid SPM & BRAINSFIT' % Hybrid SPM -> Brainsfit
        ea_docoreg_spm(moving,fixed,'nmi',0,otherfiles)
        affinefile = ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
end 
