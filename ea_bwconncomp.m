function stats=ea_bwconncomp(slicebw)
% simple wrapper for the bwconncomp issue (image processing toolbox)

try
    stats=bwconncomp(slicebw);
catch
    disp('No image processing toolbox available. Using ls_bwconncomp by Stanis?aw Adaszewski.');
    try
        stats=ls_bwconncomp(slicebw);
    catch
        try
            mex([fileparts(which('lead')),filesep,'ext_libs',filesep,'ls_conncomp_pix_list.c']);
            addpath(genpath([fileparts(which('lead')),filesep,'ext_libs',filesep]));
            stats=ls_bwconncomp(slicebw);
        catch
            error('No C-compiler found or compiling went wrong. You can solve this issue by configuring your compiler or by using a Matlab-License with image processing toolbox installed.');
            
        end
    end
end
