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
    ea_apply_normalization_tofile(options, vatspresent, wvatspresent, 1, 0);
end

if docoreg
    for vat=1:length(wvatspresent)
        copyfile(wvatspresent{vat},rwvatspresent{vat});
    end

    % BIDS FIX: Use helper function for prefix
    rest_mean = ea_prependFilename(options.prefs.rest, 'mean');
    anat_r = ea_prependFilename(options.prefs.prenii_unnormalized, 'r');
    
    reference = [directory, rest_mean];
    ea_coregimages(options, ...
    	[directory,options.prefs.prenii_unnormalized], ...
        reference, ...
        [directory, anat_r], ...
        rwvatspresent,0,[],1);
    delete([directory, anat_r]);
    ea_delete(wvatspresent);
end
