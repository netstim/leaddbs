function ea_spm_applytransforms(options,usesuit)

if ~exist('usesuit','var')
    usesuit=0;
end

directory=[options.root,options.patientname,filesep];

switch options.modality
    case 1 % MR
        postops{1}=options.prefs.prenii_unnormalized;
        postops{2}=options.prefs.tranii_unnormalized;
        postops{3}=options.prefs.cornii_unnormalized;
        postops{4}=options.prefs.sagnii_unnormalized;
        gfis{1}=options.prefs.gprenii;
        gfis{2}=options.prefs.gtranii;
        gfis{3}=options.prefs.gcornii;
        gfis{4}=options.prefs.gsagnii;
        lfis{1}=options.prefs.prenii;
        lfis{2}=options.prefs.tranii;
        lfis{3}=options.prefs.cornii;
        lfis{4}=options.prefs.sagnii;
    case 2 % CT
        postops{1}=options.prefs.prenii_unnormalized;
        postops{2}=options.prefs.ctnii_coregistered;
        gfis{1}=options.prefs.gprenii;
        gfis{2}=options.prefs.gctnii;
        lfis{1}=options.prefs.prenii;
        lfis{2}=options.prefs.ctnii;
end
[postops,gfis]=ea_appendgrid(options,postops,gfis,0);

if usesuit
postops=ea_cropsuit(postops,directory);
end
% export glfiles (a bit more coarse resolution, full brain bounding box).
for pos=1:length(gfis)
    if exist([directory,postops{pos}],'file')
        nii=ea_load_nii([directory,postops{pos}]);
        
        gaussdim=abs(nii.voxsize);
        %gaussdim=abs(gaussdim(1:3)).*2;
        if mean(gaussdim>1)
            resize_img([directory,postops{pos}],gaussdim./2,nan(2,3),0);
        else
            copyfile([directory,postops{pos}],[directory,'r',postops{pos}]);
        end
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames{1}=[directory,'r',postops{pos},''];
        matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {directory};
        matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[ea_space(options),options.primarytemplate,'.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = gaussdim;
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
        clear matlabbatch jobs;
        try movefile([directory,'swr',postops{pos}],[directory,gfis{pos}]); end
    end
end

% export lfiles (fine resolution, small bounding box).
try
    for pos=1:length(lfis)
        if exist([directory,postops{pos}],'file')
            nii=ea_load_nii([directory,postops{pos}]);
            gaussdim=abs(nii.voxsize);
            %gaussdim=abs(gaussdim(1:3)).*2;
            if mean(gaussdim>1)
                resize_img([directory,postops{pos}],gaussdim./2,nan(2,3),0);
            else
                copyfile([directory,postops{pos}],[directory,'r',postops{pos}]);
            end
            
            matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
            matlabbatch{1}.spm.util.defs.out{1}.push.fnames{1}=[directory,'r',postops{pos},''];
            matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
            matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {directory};
            matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[ea_space(options),'bb.nii']};
            matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
            matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = gaussdim;
            jobs{1}=matlabbatch;
            spm_jobman('run',jobs);
            clear matlabbatch jobs;
        end
        try movefile([directory,'swr',postops{pos}],[directory,lfis{pos}]); end
    end
end

ea_delete([directory,'rgrid.nii']);
switch options.modality
    case 1
        ea_delete([directory,'r',options.prefs.prenii_unnormalized]);
        ea_delete([directory,'r',options.prefs.tranii_unnormalized]);
        ea_delete([directory,'r',options.prefs.cornii_unnormalized]);
        ea_delete([directory,'r',options.prefs.sagnii_unnormalized]);
    case 2
        ea_delete([directory,'r',options.prefs.ctnii_coregistered]);
        ea_delete([directory,'r',options.prefs.prenii_unnormalized]);
end


function postops=ea_cropsuit(postops,directory)

[~,anatbase]=fileparts(postops{1});
% build up S for Suit:
for p=1:length(postops)
    S.resample{p}=[directory,postops{p},',1'];
end
% load affine:
load([directory,'Affine_',anatbase,'_seg1.mat']);

% load flowfield:
Vff=spm_vol([directory,'u_a_',anatbase,'_seg1.nii,1']);
% Generate grid in the space of the template
[X,Y,Z]=ndgrid(1:Vff.dim(1),1:Vff.dim(2),1:Vff.dim(3));
num_slice=Vff.dim(3);

for i=1:length(S.resample)
    [image_dir,image_name,image_ext,image_num]=spm_fileparts(S.resample{i});
    V(i)=spm_vol(S.resample{i});
    [Xm,Ym,Zm]=spmj_affine_transform(X,Y,Z,inv(Affine*V(i).mat)*Vff.mat);
    for z=1:num_slice
        Data(:,:,z)=spm_sample_vol(V(i),Xm(:,:,z),Ym(:,:,z),Zm(:,:,z),3);
    end;
    %    Data=Data.*MaskData; % Mask the image
    O(i)=V(i);
    O(i).dim=Vff.dim;
    O(i).mat=Vff.mat;
    J.images{i}{1}=fullfile(image_dir,['atemp_' image_name image_ext]);
    if (~isfield(S,'outname'))
        finalname{i}=fullfile(image_dir,['t', image_name image_ext]);
    else
        finalname{i}=S.outname{i};
    end;
    O(i).fname=J.images{i}{1};
    if (O(i).dt(1)==2)
        O(i).pinfo=[1./254*max(Data(:));0;O(i).pinfo(3)];
    end;
    spm_write_vol(O(i),Data);
    
    
    movefile([directory,postops{i}],[directory,'orig_',postops{i}]);
    movefile([directory,'atemp_',postops{i}],[directory,postops{i}]);
    
end;
