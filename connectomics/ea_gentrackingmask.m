function ea_gentrackingmask(options,threshold)
directory=[options.root,options.patientname,filesep];

ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);

%% Coreg options.prefs.prenii_unnormalized to b0 (for label.mat and FTR-Normalization)
copyfile([directory,options.prefs.prenii_unnormalized],[directory,'c',options.prefs.prenii_unnormalized]);
copyfile([directory,'c2',options.prefs.prenii_unnormalized],[directory,'cc2',options.prefs.prenii_unnormalized]);
ea_conformspaceto([directory,'c',options.prefs.prenii_unnormalized],[directory,'cc2',options.prefs.prenii_unnormalized],1,[],[],0); % make sure cc2 and anat are exactly in same space (even across software packages)
ea_backuprestore([directory,'c',options.prefs.prenii_unnormalized]);
affinefile = ea_coreg2images(options, ...
    [directory,'c',options.prefs.prenii_unnormalized], ... % moving
    [directory,options.prefs.b0], ... % fix
    [directory,'c',options.prefs.prenii_unnormalized], ... % out
    {[directory,'cc2',options.prefs.prenii_unnormalized]}, ... % other
    1); % writeout transform

% Rename saved transform for further use: canat2b0*.mat to anat2b0*.mat,
% and b02canat2b0*.mat to b02anat*.mat
[~, anat] = ea_niifileparts(options.prefs.prenii_unnormalized);
if ~isempty(affinefile)
    for i = 1:numel(affinefile)
        movefile(affinefile{i}, strrep(affinefile{i}, ['c',anat], anat));
    end
end

movefile([directory,'cc2',options.prefs.prenii_unnormalized],[directory,'trackingmask.nii']);
delete([directory,'c',options.prefs.prenii_unnormalized]);

tr=ea_load_nii([options.root,options.patientname,filesep,'trackingmask.nii']);
if threshold
    tr.img=tr.img>0.1;
    tr.fname=[options.root,options.patientname,filesep,'ttrackingmask.nii'];
    ea_write_nii(tr);
end
