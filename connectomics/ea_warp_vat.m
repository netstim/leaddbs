function ea_warp_vat(ref_filename,b0rest,options,handles)
directory=[options.root,options.patientname,filesep];

stims=get(handles.vatseed,'String');
stim=stims{get(handles.vatseed,'Value')};

% check which vat-files are present:
vatfnames={[directory,'stimulations',filesep,stim,filesep,'','vat_right.nii']
    [directory,'stimulations',filesep,stim,filesep,'','vat_left.nii']};
cnt=1;
donorm=0;
docoreg=0;
for vatfname=1:2
    if exist(vatfnames{vatfname},'file')
        vatspresent{cnt}=vatfnames{vatfname};
        [pth,fn,ext]=fileparts(vatfnames{vatfname});
        wvatspresent{cnt}=[pth,filesep,'w',fn,ext];
        if ~exist(wvatspresent{cnt},'file')
            donorm=1;
        end
        rwvatspresent{cnt}=[pth,filesep,'r',b0rest,'w',fn,ext];
        if ~exist(rwvatspresent{cnt},'file')
            docoreg=1;
        end
        cnt=cnt+1;
    end
end


    
if donorm
    %% warp vat into pre_tra-space:  
    

           
    switch spm('ver')
        case 'SPM8'
            
    
            matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']};
            matlabbatch{1}.spm.util.defs.ofname = '';
            matlabbatch{1}.spm.util.defs.fnames = vatspresent;
            matlabbatch{1}.spm.util.defs.savedir.saveusr = {[directory,'stimulations',filesep,stim,filesep,filesep]};
            matlabbatch{1}.spm.util.defs.interp = 0;
     
        case 'SPM12'
            
            matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']};
            matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = vatspresent;
            matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {[directory,'stimulations',filesep,stim,filesep,filesep]};
            matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
            matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            
    end
    % execute batch..
    cfg_util('run',{matlabbatch});
    clear matlabbatch
end
if docoreg
    %% coreg vat into b0/rest-space:
    copyfile([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized]);
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[options.root,options.patientname,filesep,ref_filename,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = wvatspresent;
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = ['r',b0rest];
    cfg_util('run',{matlabbatch});
    clear matlabbatch
    
    delete([options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized]);
    delete([options.root,options.patientname,filesep,'r',b0rest,'c',options.prefs.prenii_unnormalized]);
end