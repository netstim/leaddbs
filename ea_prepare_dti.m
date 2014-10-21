function ea_prepare_dti(options)
% calculates diffusion tensor etc. using the Freiburg DTI&Fibertools for
% SPM8.

    disp('Building DTI files...');
   load([options.root,options.patientname,filesep,options.prefs.bval]);
   [~,bvfname]=fileparts(options.prefs.bval);
   bvals=eval(bvfname);
   ea_build_DTD(max(bvals),[options.root,options.patientname,filesep],options.prefs.dti,options.prefs.DTD,options.prefs.HARDI,options.prefs.bval,options.prefs.bvec);

   % export B0
   matlabbatch{1}.impexp_NiftiMrStruct.bo2nifti.srcdtdchoice.srcdtdstruct = {[options.root,options.patientname,filesep,options.prefs.DTD]};
   matlabbatch{1}.impexp_NiftiMrStruct.bo2nifti.outname.outimg.outdir = {[options.root,options.patientname,filesep]};
   matlabbatch{1}.impexp_NiftiMrStruct.bo2nifti.outname.outimg.fname = options.prefs.b0;
   jobs{1}=matlabbatch;
   cfg_util('run',jobs);
   clear jobs matlabbatch
   disp('Done.');