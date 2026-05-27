function [AllX,space]=ea_unifiedmapping_exportefieldmap(vatlist,obj)

disp('Need to export Efields in proper format, this may take a while');
prefs=ea_prefs;
outdir = fullfile(fileparts(obj.leadgroup), 'sweetspots', obj.ID, filesep);
ea_mkdir(outdir);
if size(vatlist,2)>1
    sidesuffices = {'_r', '_l'};
else
    sidesuffices={''};
end

% create common space across niftis
space = cell(1, size(vatlist,2));
for side=1:size(vatlist,2)
    fname = [outdir, 'efield_bb', sidesuffices{side}, '.nii'];
    if ~isfield(obj.M,'pseudoM')
        ea_create_threshold_union_ref_nii(fname, vatlist(isfile(vatlist(:,side)),side), prefs.explorer.efieldlowfilter, repmat(obj.statsettings.sweetspotresolution,1,3)); % at least one efield should have a value of 150 for this to be included
    else
        ea_create_union_ref_nii([outdir, 'efield_bb.nii'], vatlist(isfile(vatlist(:,side)),side), obj.roithresh, repmat(obj.statsettings.sweetspotresolution,1,3));
        ea_cprintf('CmdWinWarnings', "PseudoM VTAs are used, low field filtering disabled! This can result in a large RAM consumption.\n");
    end
    nii=ea_load_nii(fname); % reload for space function.
    space{side}=nii;
end

% now conform each VTA to space
AllX=cell(size(vatlist,2),1);
for vat=1:size(vatlist,1)
    for side=1:size(vatlist,2)
        if isfile(vatlist{vat,side})
            id=ea_generate_uuid;
    
            copyfile(vatlist{vat,side},[outdir,'tmp_efield',id,'.nii']);
    
            ea_fsl_reslice([outdir,'tmp_efield',id,'.nii'],[outdir,'efield_bb',sidesuffices{side},'.nii'],[outdir,'tmp_efield',id,'.nii'],'nearestneighbour');
    
            nii=ea_load_nii([outdir,'tmp_efield',id,'.nii']);
    
            if ~exist('AllX','var')
                AllX{side}=nan(size(vatlist,1),numel(nii.img));
            end
            AllX{side}(vat,:)=nii.img(:);
            ea_delete([outdir,'tmp_efield',id,'.nii']);
        end
    end
end

ea_delete([outdir,'efield_bb.nii']);
ea_delete([outdir,'efield_bb_l.nii']);
ea_delete([outdir,'efield_bb_r.nii']);
ea_delete([outdir,'bb_nan.nii']);
ea_delete(outdir);
