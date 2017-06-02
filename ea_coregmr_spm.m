function ea_coregmr_spm(options,doreslice)
% this function coregisters postoperative images to preoperative images
% using SPM.


directory = [options.root,options.prefs.patientdir];

if ~exist([directory,filesep,options.prefs.prenii_unnormalized],'file')
    warning('No preoperative acquisition found. Coregistration not possible.');
    return
end

costfuns={'nmi','mi','ecc','ncc'};
cfundo=[2,1];
for export=1:3
    for costfun=1:length(cfundo)
        switch export
            case 1
                fina=[directory,filesep,options.prefs.tranii_unnormalized,',1'];
            case 2
                fina=[directory,filesep,options.prefs.cornii_unnormalized,',1'];
            case 3
                fina=[directory,filesep,options.prefs.sagnii_unnormalized,',1'];
        end
        if exist(fina(1:end-2),'file')
            if costfun==length(cfundo) % only at final stage apply refining if set
                
                ea_docoreg_spm(options,fina,[directory,filesep,options.prefs.prenii_unnormalized,',1'], ...
                               costfuns{cfundo(costfun)},doreslice,{''},options.prefs.mrcoreg.writeoutcoreg);
            else % dont reslice, dont refine (not last pass).
                ea_docoreg_spm(options,fina,[directory,filesep,options.prefs.prenii_unnormalized,',1'],...
                               costfuns{cfundo(costfun)},0,{''},options.prefs.mrcoreg.writeoutcoreg);
            end
            disp(['*** Coregistration pass (',costfuns{cfundo(costfun)},') completed.']);
        end
        
    end
end

if doreslice
    try   movefile([directory,filesep,'r',options.prefs.tranii_unnormalized],[directory,filesep,options.prefs.tranii_unnormalized]); end
    try   movefile([directory,filesep,'r',options.prefs.cornii_unnormalized],[directory,filesep,options.prefs.cornii_unnormalized]); end
    try   movefile([directory,filesep,'r',options.prefs.sagnii_unnormalized],[directory,filesep,options.prefs.sagnii_unnormalized]); end
end
