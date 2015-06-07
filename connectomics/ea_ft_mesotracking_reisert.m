function varargout=ea_ft_mesotracking_reisert(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Mesoscopic Fibertracking (Reisert et al. 2014).';
    varargout{2}={'SPM8','SPM12'};
    return
end

gdti_trackingparams='standard'; % select param-preset here (see below, also to create your own)


switch gdti_trackingparams
    
    case 'hd_book'
        para_weight = 0.0006;
        para_other = [1
            0.001
            50
            5000000000
            0.5
            3.75
            0.2
            1];
    case 'hd_a'
        para_weight = 0.006;
        para_other = [1
            0.001
            50
            5000000000
            0.5
            3.75
            0.2
            1];
    case 'hd_a_light'
        para_weight = 0.02;
        para_other = [1
            0.001
            50
            300000000
            0.5
            3
            0.2
            1];
    case 'hd_a_verylight'
        para_weight = 0.03;
        para_other = [0.1
            0.001
            50
            300000000
            0.5
            3
            0.2
            1];
    case 'standard'
        para_weight = 0.058;
        para_other = [0.1
            0.001
            50
            300000000
            1
            3
            0.2
            1];
    case 'standard_enhanced'
        para_weight = 0.058;
        para_other = [0.1
            0.001
            50
            300000000
            1
            3
            0.2
            1.5];
        
end

keyboard
directory=[options.root,options.patientname,filesep];
ea_exportb0(options);

% create c2 from anat
ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
%% coreg anat to b0

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,options.prefs.b0,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,'c2',options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
cfg_util('run',{matlabbatch}); clear matlabbatch



%% mesoft part goes here

mesoGT_tool('loadData','nii',[directory,options.prefs.dti],{[directory,options.prefs.bval],[directory,options.prefs.bval]},[directory,'rc2',options.prefs.prenii_unnormalized],0.5);
% 
% 
% %% export .trk copy for trackvis visualization
% 
% dnii=load_nii([directory,options.prefs.b0]);
% niisize=size(dnii.img); % get dimensions of reference template.
% specs.origin=[0,0,0];
% specs.dim=niisize;
% try
%     H=spm_dicom_headers([directory,prefs.sampledtidicom]);
%     specs.orientation=H{1,1}.ImageOrientationPatient;
% catch
%     specs.orientation=[1,0,0,0,1,0];
% end
% [~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
% ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
% disp('Done.');


function ea_exportb0(options)

bvals=load([options.root,options.patientname,filesep,prefs.bval]);
idx=find(bvals==0);

for fi=idx
   fis{fi}=[options.root,options.patientname,filesep,prefs.dti,',',num2str(fi)];
end

matlabbatch{1}.spm.util.imcalc.input = fis;
matlabbatch{1}.spm.util.imcalc.output = [options.prefs.b0];
matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname]};
matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
cfg_util('run',{matlabbatch}); clear matlabbatch


