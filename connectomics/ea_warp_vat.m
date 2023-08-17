function ea_warp_vat(b0rest, options, handles)
directory=[options.root,options.patientname,filesep];

if strcmp(b0rest,'rest') % processing rest files
    b0rest=ea_stripext(options.prefs.rest);
end

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

        rwvatspresent{cnt}=[pth,filesep,'r',b0rest,'w',fn,'.nii'];
        if ~exist(rwvatspresent{cnt},'file')
            docoreg=1;
        end

        cnt=cnt+1;
    end
end

if donorm
    %% warp vat into pre_tra-space:
    whichnormmethod=ea_whichnormmethod([directory]);
    switch whichnormmethod
        case ea_getantsnormfuns

            ea_ants_apply_transforms(options, ...
                vatspresent, ...
                wvatspresent,...
                1,'','','NearestNeighbor');

        case ea_getfslnormfuns

            ea_fsl_apply_normalization(options, ...
                vatspresent, ...
                wvatspresent,...
                1,'','','nn');
        otherwise

            matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
            matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = vatspresent;
            matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {[directory,'stimulations',filesep,stim,filesep,filesep]};
            matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
            matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];

            % execute batch
            spm_jobman('run',{matlabbatch});
            clear matlabbatch
    end
end

if docoreg
    for vat=1:length(wvatspresent)
        copyfile(wvatspresent{vat},rwvatspresent{vat});
    end

    reference = [directory,'mean',options.prefs.rest];
    ea_coregimages(options, ...
    	[directory,options.prefs.prenii_unnormalized], ...
        reference, ...
        [directory,'r',options.prefs.prenii_unnormalized], ...
        rwvatspresent,0,[],1);
    delete([directory,'r',options.prefs.prenii_unnormalized]);
    ea_delete(wvatspresent);
end
