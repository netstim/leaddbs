function ea_coreg2images(options,moving,fixed,ofile)



[directory,mfilen,ext]=fileparts(moving);
directory=[directory,filesep];
mfilen=[mfilen,ext];

if exist([directory,'raw_',mfilen],'file');
    copyfile([directory,'raw_',mfilen],[directory,mfilen]);
end
copyfile([directory,mfilen],[directory,'raw_',mfilen]);

switch options.coregmr.method
    case 1 % SPM
        ea_docoreg_spm(moving,fixed,'nmi',1)
        movefile([directory,'r',mfilen],ofile);
    case 2 % ANTs
        ea_ants(fixed,...
            moving,...
            ofile,0);
    case 3 % BRAINSFit
        ea_brainsfit(fixed,...
            moving,...
            ofile,0);
    case 4 % Hybrid SPM -> ANTs
        ea_docoreg_spm(moving,fixed,'nmi',0)
        ea_ants(fixed,...
            moving,...
            ofile,0);
    case 5 % Hybrid SPM -> Brainsfit
        ea_docoreg_spm(moving,fixed,'nmi',0)
        ea_brainsfit(fixed,...
            moving,...
            ofile,0);
end