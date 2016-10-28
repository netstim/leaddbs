function ea_coreg2images(options,moving,fixed,ofile,otherfiles,writeoutmat)

if ~exist('otherfiles','var')
   otherfiles={''}; 
end
if ~exist('writeoutmat','var')
    writeoutmat=0;
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
        
        ea_docoreg_spm(moving,fixed,'nmi',1,otherfiles)
        movefile([directory,'r',mfilen],ofile);
        for ofi=1:length(otherfiles)
           [pth,fn,ext]=fileparts(otherfiles{ofi});
           try % could be empty cell.
           movefile(fullfile(pth,['r',fn,ext]),fullfile(pth,[fn,ext]));
           end
        end
        
    case 'Coreg MRIs: ANTs' % ANTs
        
        ea_ants(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Coreg MRIs: BRAINSFIT' % BRAINSFit
        ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Coreg MRIs: Hybrid SPM & ANTs' % Hybrid SPM -> ANTs
        ea_docoreg_spm(moving,fixed,'nmi',0,otherfiles)
        ea_ants(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
    case 'Coreg MRIs: Hybrid SPM & BRAINSFIT' % Hybrid SPM -> Brainsfit
        ea_docoreg_spm(moving,fixed,'nmi',0,otherfiles)
        ea_brainsfit(fixed,...
            moving,...
            ofile,writeoutmat,otherfiles);
end