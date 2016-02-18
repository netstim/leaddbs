function [finas]=ea_coregmr(options,automan)

if options.modality==2 % in CT imaging, coregistration is done elsewhere.
    return
end

if strcmp(options.prefs.mrcoreg.default,'ants')
   ea_coregmr_ants(options);
   return 
end



if exist([options.root,options.patientname,filesep,'raw_',options.prefs.tranii_unnormalized],'file')
    copyfile([options.root,options.patientname,filesep,'raw_',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.tranii_unnormalized]);
end
if exist([options.root,options.patientname,filesep,'raw_',options.prefs.cornii_unnormalized],'file')
    copyfile([options.root,options.patientname,filesep,'raw_',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.cornii_unnormalized]);
end
if exist([options.root,options.patientname,filesep,'raw_',options.prefs.sagnii_unnormalized],'file')
    copyfile([options.root,options.patientname,filesep,'raw_',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized]);
end


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


for export=1:2
    for costfun=cfundo
        switch export
            case 1
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized,',1'];
            case 2
                fina=[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized,',1'];
        end


        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {fina};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = costfuns{costfun};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [12 10 8 6 4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

        jobs{1}=matlabbatch;
        try

            cfg_util('run',jobs);
            normlog(1)=1;
            disp('*** Coregistration between transversal and coronar versions worked.');
            finas{export}=fina; % assign only if worked.
        catch
            disp('*** Coregistration between transversal and coronar versions failed / Using CT Modality.');
            %ea_error('This normalization cannot be performed automatically with eAuto. Try using different software for the normalization step. Examples are to use SPM directly, or to use FSL, Slicer or Bioimaging Suite.');
        end
        clear matlabbatch jobs;
        if manual
            matlabbatch{1}.spm.util.checkreg.data = {[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1'];
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




% second step, coreg post to pre.

normlog=zeros(4,1); % log success of processing steps. 4 steps: 1. coreg tra and cor, 2. grand mean normalization 3. subcortical normalization 4. subcortical fine normalization that spares the ventricles.



for costfun=cfundo


    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']};
    try
        matlabbatch{1}.spm.spatial.coreg.estimate.other = finas;
    catch
        matlabbatch{1}.spm.spatial.coreg.estimate.other={''};
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = costfuns{costfun};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [15 10 8 6 4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [8 8];

    jobs{1}=matlabbatch;
    try
        cfg_util('run',jobs);
        normlog(1)=1;
        disp('*** Coregistration between pre and post versions worked.');
    catch
        disp('*** Coregistration between pre and post versions failed.');
        %ea_error('This normalization cannot be performed automatically with LEAD. Try using different software for the normalization step. Examples are to use SPM directly, or to use FSL, Slicer or Bioimaging Suite.');
    end
    clear matlabbatch jobs;

    if manual
        if exist('finas','var')
            matlabbatch{1}.spm.util.checkreg.data = [{[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']},finas];
        else
            matlabbatch{1}.spm.util.checkreg.data = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1'],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,',1']};

        end
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear matlabbatch jobs;


        yninp = input('Please check reg between Pre- and Post-OP versions. Is result precise? (y/n)..','s');
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
end

if ~exist('finas','var')
    finas={};
end
finas=[finas,{[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']}];


