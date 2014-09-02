function convertfromTPMtoicbm2009(folder)
% enter the root folder of your patients directory and normalized files as
% well as reconstructions will be converted from SPM's TPM.nii space to the
% ICBM 2009 b template space (which is used by LEAD and can be found in
% lead_dbs/templates/mni_hires.nii
% Please note: The current release of LEAD-DBS uses a TPM.nii file that
% already *is* in the new ICBM space. This script should thus only be used
% for reconstructions made with older releases of LEAD-DBS and only with
% normalizations performed with the "SPM New Segment" approach. All other
% normalization approaches have normalized to the ICBM 2009b release ever
% since.

if ~strcmp(folder(end),filesep)
    folder=[folder,filesep];
end



tmat=[    1.0380   -0.0055    0.0021   -0.9238
    0.0080    1.0394   -0.0045    1.0272
   -0.0057   -0.0101    1.0162    0.1113
         0         0         0    1.0000];


load([folder,'ea_reconstruction.mat'])

% convert coords_mm

for side=1:2
    tc=[coords_mm{side},ones(size(coords_mm{side},1),1)];
    coords_mm{side}=tmat*tc';
    coords_mm{side}=coords_mm{side}(1:3,:)';
end

% convert trajectory

for side=1:2
    tc=[trajectory{side},ones(size(trajectory{side},1),1)];
    trajectory{side}=tmat*tc';
    trajectory{side}=trajectory{side}(1:3,:)';
end

save([folder,'ea_reconstruction.mat'],'trajectory','coords_mm');


% transform images:

images={'ltra.nii','lcor.nii','lsag.nii','gltra.nii','glcor.nii','glsag.nii','lpre.nii','glpre.nii','lfusion.nii','glfusion.nii'};

for image=1:length(images)
    if ~exist([folder,images{image}],'file') && exist([folder,images{image},'.gz'],'file')
        gunzip([folder,images{image},'.gz'])
        wasgz=1;
    else
        wasgz=0;
    end
    
    matlabbatch{1}.spm.util.reorient.srcfiles = {[folder,images{image},',1']};
    matlabbatch{1}.spm.util.reorient.transform.transM = [1.03799350893443 -0.00551554932118825 0.00211932303236285 -0.923782191169304
        0.00800793946529943 1.0393803580894 -0.00445193813236372 1.02715996003749
        -0.00568315298822271 -0.0101344712041121 1.01619757821906 0.111266550723032
        0 0 0 1];
    matlabbatch{1}.spm.util.reorient.prefix = '';
    jobs{1}=matlabbatch;
    try
        cfg_util('run',jobs);
    end
    clear jobs matlabbatch
    
    if wasgz
       gzip([folder,images{image}]) 
    end
end

