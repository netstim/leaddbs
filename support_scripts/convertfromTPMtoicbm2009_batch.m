function convertfromTPMtoicbm2009_batch(ptroot)
% enter a patient folder and normalized files as
% well as reconstructions will be converted from SPM's TPM.nii space to the
% ICBM 2009 b template space (which is used by LEAD and can be found in
% lead_dbs/templates/mni_hires.nii
% Please note: The current release of LEAD-DBS uses a TPM.nii file that
% already *is* in the new ICBM space. This script should thus only be used
% for reconstructions made with older releases of LEAD-DBS and only with
% normalizations performed with the "SPM New Segment" approach. All other
% normalization approaches have normalized to the ICBM 2009b release ever
% since.

root=[fileparts(ptroot),filesep];
pts=dir(ptroot);

for pt=1:length(pts)
    if pts(pt).isdir
       try
           load([root,pts(pt).name,filesep,'ea_normmethod_applied.mat']);
           if strcmp(norm_method_applied,'ea_normalize_spmnewseg');
               convertfromTPMtoicbm2009([root,pts(pt).name,filesep]);
           end
       end 
    end
end 