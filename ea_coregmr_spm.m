function ea_coregmr_spm(options,doreslice,refine)
% this function coregisters postoperative images to preoperative images
% using SPM.

costfuns={'nmi','mi','ecc','ncc'};
cfundo=[2,1];
for export=1:3
    cnt=1;
    for costfun=cfundo
        switch export
            case 1
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1'];
            case 2
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized,',1'];
            case 3
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized,',1'];
        end
        if exist(fina(1:end-2),'file')
            if cnt==length(cfundo) % only at final stage apply refining if set
                ea_docoreg_spm(options,fina,[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'],costfuns{costfun},doreslice,{''},0,refine);
            else % dont reslice, dont refine (not last pass).
                ea_docoreg_spm(options,fina,[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'],costfuns{costfun},0,{''},0,0);
            end
            cnt=cnt+1;
            disp(['*** Coregistration pass (',costfuns{costfun},') completed.']);
        end
        
    end
end

if doreslice
    try   movefile([options.root,options.patientname,filesep,'r',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.tranii_unnormalized]); end
    try   movefile([options.root,options.patientname,filesep,'r',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.cornii_unnormalized]); end
    try   movefile([options.root,options.patientname,filesep,'r',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized]); end
end
