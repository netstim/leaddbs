function locked = ea_coreglocked(options, moving)
locked=0;

if isfield(options, 'overwriteapproved') && options.overwriteapproved
    return
end

directory=[options.root, options.patientname, filesep];
if ~exist([directory,'ea_coreg_approved.mat'], 'file')
    return
else
    approved = load([directory, 'ea_coreg_approved.mat']);
    [~, moving] = ea_niifileparts(moving);
    if isfield(approved, moving)
        locked = approved.(moving);
    end
end
