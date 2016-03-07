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
            ea_docoreg(fina,[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'],costfuns{costfun},doreslice)
            normlog(1)=1;
            disp('*** Coregistration between transversal and coronar versions worked.');
            finas{export}=fina; % assign only if worked.
        catch
            disp('*** Coregistration between transversal and coronar versions failed / Using CT Modality.');
            %ea_error('This normalization cannot be performed automatically with eAuto. Try using different software for the normalization step. Examples are to use SPM directly, or to use FSL, Slicer or Bioimaging Suite.');
        end
        if manual
            matlabbatch{1}.spm.util.checkreg.data = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'];
                fina};
            jobs{1}=matlabbatch;
            try % CT
                cfg_util('run',jobs);

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


function ea_docoreg(moving,fixed,cfun,doreslice)
if doreslice
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fixed};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {moving};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = cfun;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
else
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fixed};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {moving};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = cfun;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [12 10 8 6 4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
end

cfg_util('run',{matlabbatch});