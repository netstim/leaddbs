function ea_prepare_dti(options)
% calculates diffusion tensor etc. using the Freiburg DTI&Fibertools for
% SPM8. If DTI&Fibertools are not installed on your computer, please obtain
% them from here:
% http://www.uniklinik-freiburg.de/mr-en/research-groups/diffperf/fibertools.html
% and add them to your Matlab-path to install them.

ea_checkdti_ft % checks if dti_tool is available (DTI&Fibertools).

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



function ea_checkdti_ft
if ~exist('dti_tool','file')==2
    ea_error('Please install the DTI&Fibertools for SPM to run DTI-Fibertracking (http://www.uniklinik-freiburg.de/mr-en/research-groups/diffperf/fibertools.html).')
end
