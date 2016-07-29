function ea_coregmr_spm(options,automan,doreslice)

costfuns={'nmi','mi','ecc','ncc'};

% first step, coregistration between transversal and coronar/sagittal versions. on full brain

switch automan
    case 'manual'
        cfundo=1:4;
        manual=1;
    case 'auto'
        cfundo=[2,1,3,4];
        manual=0;
    otherwise
        ea_error('Coregistration prefs must be either set to auto or manual. Please modify ea_prefs.m accordingly.');
end

normlog=zeros(4,1); % log success of processing steps. 4 steps: 1. coreg tra and cor, 2. grand mean normalization 3. subcortical normalization 4. subcortical fine normalization that spares the ventricles.

for export=1:3
    for costfun=cfundo
        switch export
            case 1
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1'];
            case 2
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized,',1'];
            case 3
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized,',1'];
        end

        try
            if exist(fina(1:end-2),'file')
                ea_docoreg_spm(fina,[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'],costfuns{costfun},doreslice)
                normlog(1)=1;
                disp(['*** Coregistration between transversal and coronar versions worked (',costfuns{costfun},').']);
                finas{export}=fina; % assign only if worked.
            end
        catch
            disp('*** Coregistration between transversal and coronar versions failed / Using CT Modality.');
            %ea_error('This normalization cannot be performed automatically with eAuto. Try using different software for the normalization step. Examples are to use SPM directly, or to use FSL, Slicer or Bioimaging Suite.');
        end
        
        if manual
            matlabbatch{1}.spm.util.checkreg.data = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'];
                fina};
            jobs{1}=matlabbatch;
            try % CT
                spm_jobman('run',jobs);

                yninp = input('Please check reg between Post-OP versions. Is result precise? (y/n)..','s');
                if strcmpi(yninp,'y')
                    disp('Good. Moving on...');
                    break
                else
                    if costfun==4
                        ea_error('Problem cannot be solved automatically.')
                    else
                        disp('Trying with another cost-function');
                    end
                end
            end
            clear matlabbatch jobs;
        end
    end
end

if doreslice
    try   movefile([options.root,options.patientname,filesep,'r',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.tranii_unnormalized]); end
    try   movefile([options.root,options.patientname,filesep,'r',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.cornii_unnormalized]); end
    try   movefile([options.root,options.patientname,filesep,'r',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized]); end
end
