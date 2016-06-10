function ea_prepare_dti(options)
% calculates diffusion tensor etc. using the Freiburg DTI&Fibertools
% http://www.uniklinik-freiburg.de/mr-en/research-groups/diffperf/fibertools.html

if ~exist([options.root,options.patientname,filesep,options.prefs.HARDI],'file');
    disp('Building DTI files...');
    
%    load([options.root,options.patientname,filesep,options.prefs.bval]);
%    [~,bvfname]=fileparts(options.prefs.bval);
%    bvals=eval(bvfname);
%    ea_build_DTD(max(bvals),[options.root,options.patientname,filesep],options.prefs.dti,options.prefs.DTD,options.prefs.HARDI,options.prefs.bval,options.prefs.bvec);
%   

try %unring
    dti=ea_load_untouch_nii([options.root,options.patientname,filesep,options.prefs.dti]);
    dti.img=ea_unring(dti.img);
    ea_save_untouch_nii(dti,[options.root,options.patientname,filesep,options.prefs.dti]);
end
    
    % build HARDI in new way:
    
    
    hr=ea_DPS_nifti_to_hardi([options.root,options.patientname,filesep,options.prefs.dti],...
        [options.root,options.patientname,filesep,options.prefs.bvec],...
        [options.root,options.patientname,filesep,options.prefs.bval]);
    
    mrstruct_write(hr,[options.root,options.patientname,filesep,options.prefs.HARDI])
 
    
%     dtd=ea_convertHARDI2DTD(dtdstruct_read([options.root,options.patientname,filesep,options.prefs.HARDI]));
%     dtstruct_write(dtd,[options.root,options.patientname,filesep,options.prefs.DTD])
%     
    % export B0
    if ~exist([options.root,options.patientname,filesep,options.prefs.b0],'file');
        ea_exportb0(options);
    end

else
    disp('HARDI found, no need to rebuild.');
end