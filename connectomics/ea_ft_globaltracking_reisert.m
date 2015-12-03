function varargout=ea_ft_globaltracking_reisert(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Global Fibertracking (Reisert et al. 2011)';
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



%% build DTD (tensor calculation)


ea_prepare_dti(options)

%% mask for tracking

directory=[options.root,options.patientname,filesep];


% 'new segment' options.prefs.prenii_unnormalized
ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);


%% Coreg options.prefs.prenii_unnormalized to b0 (for label.mat and FTR-Normalization)
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,options.prefs.b0,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[directory,'c1',options.prefs.prenii_unnormalized,',1'];
    [directory,'c2',options.prefs.prenii_unnormalized,',1'];
    [directory,'c3',options.prefs.prenii_unnormalized,',1'];
    [directory,'c4',options.prefs.prenii_unnormalized,',1'];
    [directory,'c5',options.prefs.prenii_unnormalized,',1']
    };
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'rb0';

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;

%% build tracking mask

matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.srcimgs = {
    [directory,'rb0c1',options.prefs.prenii_unnormalized,',1'];
    [directory,'rb0c2',options.prefs.prenii_unnormalized,',1'];
    [directory,'rb0c3',options.prefs.prenii_unnormalized,',1'];
    [directory,'rb0c4',options.prefs.prenii_unnormalized,',1'];
    [directory,'rb0c5',options.prefs.prenii_unnormalized,',1']
    };
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.thresh = [0.5 Inf];
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.nfval = false;
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.outchoice.outmat.outdir = {directory};
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.outchoice.outmat.fname = 'c1_5_mask';


jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;


load([directory,'c1_5_mask.mat']);
delete([directory,'c1_5_mask.mat']);

[~,pretra]=fileparts(options.prefs.prenii_unnormalized);

%invert C3-5
maskCell{length(maskCell)+1}=~maskCell{find(ismember(maskNamesCell,['rb0c3',pretra])==1)};
maskNamesCell{length(maskNamesCell)+1}=['!rb0c3',pretra];
maskCell{length(maskCell)+1}=~maskCell{find(ismember(maskNamesCell,['rb0c4',pretra])==1)};
maskNamesCell{length(maskNamesCell)+1}=['!rb0c4',pretra];
maskCell{length(maskCell)+1}=~maskCell{find(ismember(maskNamesCell,['rb0c5',pretra])==1)};
maskNamesCell{length(maskNamesCell)+1}=['!rb0c5',pretra];

%add !c3-5 to c2 (noise reduction)
maskCell{length(maskCell)+1}=maskCell{find(ismember(maskNamesCell,['!rb0c3',pretra])==1)}.*maskCell{find(ismember(maskNamesCell,['rb0c2',pretra])==1)};
maskNamesCell{length(maskNamesCell)+1}=['c2*!rb0c3',pretra];
maskCell{length(maskCell)+1}=maskCell{find(ismember(maskNamesCell,['!rb0c4',pretra])==1)}.*maskCell{find(ismember(maskNamesCell,['c2*!rb0c3',pretra])==1)};
maskNamesCell{length(maskNamesCell)+1}=['c2*!c3*!rb0c4',pretra];
maskCell{length(maskCell)+1}=maskCell{find(ismember(maskNamesCell,['!rb0c5',pretra])==1)}.*maskCell{find(ismember(maskNamesCell,['c2*!c3*!rb0c4',pretra])==1)};
maskNamesCell{length(maskNamesCell)+1}='trackingmask';

maskCell{length(maskCell)+1}=maskCell{find(ismember(maskNamesCell,['rb0c1',pretra])==1)}+maskCell{find(ismember(maskNamesCell,['trackingmask'])==1)};
maskNamesCell{length(maskNamesCell)+1}='greymask';

trackingmask.maskCell=maskCell;
trackingmask.maskNamesCell=maskNamesCell;
trackingmask.mrsProp=mrsProp;
trackingmask.sizeAy=sizeAy;
trackingmask.user=user;
trackingmask.version=version;

trackingmask=maskstruct_modify(trackingmask,'dilatation','trackingmask','diltrackingmask','XYZ',5);



save([directory,'trackingmask'],'-struct','trackingmask');


%% do tracking



matlabbatch{1}.dtijobs.tracking.GTtrack.fname.filenameHARDI = {[directory,options.prefs.HARDI]};
matlabbatch{1}.dtijobs.tracking.GTtrack.newfile.out.dir = {directory};
matlabbatch{1}.dtijobs.tracking.GTtrack.newfile.out.fname = options.prefs.FTR_unnormalized;
matlabbatch{1}.dtijobs.tracking.GTtrack.trackingarea.maskstruct.filenameMASK = {[directory,'trackingmask.mat']};
matlabbatch{1}.dtijobs.tracking.GTtrack.trackingarea.maskstruct.roiid = 'trackingmask';
matlabbatch{1}.dtijobs.tracking.GTtrack.parameters = 1;
matlabbatch{1}.dtijobs.tracking.GTtrack.para_weight.custom_para_weight = para_weight;
matlabbatch{1}.dtijobs.tracking.GTtrack.para_other.custom_para_other = para_other;
matlabbatch{1}.dtijobs.tracking.GTtrack.minlen = 3;
matlabbatch{1}.dtijobs.tracking.GTtrack.maxlen = Inf;


jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;





%% export .trk copy for trackvis visualization

dnii=load_nii([directory,options.prefs.b0]);
niisize=size(dnii.img); % get dimensions of reference template.
specs.origin=[0,0,0];
specs.dim=niisize;
try
    H=spm_dicom_headers([directory,prefs.sampledtidicom]);
    specs.orientation=H{1,1}.ImageOrientationPatient;
catch
    specs.orientation=[1,0,0,0,1,0];
end
[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
disp('Done.');





