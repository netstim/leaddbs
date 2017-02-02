function ea_coregmr_spm(options,doreslice,refine)
% this function coregisters postoperative images to preoperative images
% using SPM.

directory = [options.root,options.prefs.patientdir];
costfuns={'nmi','mi','ecc','ncc'};
cfundo=[2,1];
for export=1:3
    cnt=1;
    for costfun=cfundo
        switch export
            case 1
                fina=[directory,filesep,options.prefs.tranii_unnormalized,',1'];
            case 2
                fina=[directory,filesep,options.prefs.cornii_unnormalized,',1'];
            case 3
                fina=[directory,filesep,options.prefs.sagnii_unnormalized,',1'];
        end
        if exist(fina(1:end-2),'file')
            if cnt==length(cfundo) % only at final stage apply refining if set
                ea_docoreg_spm(options,fina,[directory,filesep,options.prefs.prenii_unnormalized,',1'], ...
                               costfuns{costfun},doreslice,{''},options.prefs.mrcoreg.writeoutcoreg,refine);
            else % dont reslice, dont refine (not last pass).
                ea_docoreg_spm(options,fina,[directory,filesep,options.prefs.prenii_unnormalized,',1'],...
                               costfuns{costfun},0,{''},options.prefs.mrcoreg.writeoutcoreg,0);
            end
            cnt=cnt+1;
            disp(['*** Coregistration pass (',costfuns{costfun},') completed.']);
        end
        
    end
end

if doreslice
    try   movefile([directory,filesep,'r',options.prefs.tranii_unnormalized],[directory,filesep,options.prefs.tranii_unnormalized]); end
    try   movefile([directory,filesep,'r',options.prefs.cornii_unnormalized],[directory,filesep,options.prefs.cornii_unnormalized]); end
    try   movefile([directory,filesep,'r',options.prefs.sagnii_unnormalized],[directory,filesep,options.prefs.sagnii_unnormalized]); end
end
