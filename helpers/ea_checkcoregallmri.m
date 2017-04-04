function ea_checkcoregallmri(options,usebrainmask,usefa)    
if ~exist('usebrainmask','var')
    usebrainmask=0; % by default dont use brain mask.
end
if ~exist('usefa','var')
    usefa=1; % by default use it.
end

if ~ea_seemscoregistered(options) % check headers of files to see if already coregistered.
    if strcmp(options.coregmr.method,'Do not coregister MRIs (already coregistered)')
        ea_error('You specified that images are already coregistered & resliced but based on their NIFtI Headers they seem not. Please set an MR coregistration option and re-run.');
    end
    ea_coreg_all_mri(options,usebrainmask)
end
if usefa
    ea_coreg_fa(options,usebrainmask);
end