function ea_compat_pt(uipatdir)

for pt=1:length(uipatdir)
    directory=[uipatdir{pt},filesep];
    
    fids=dir([directory,'fiducials',filesep,'*.nii.gz']);
    for fid=1:length(fids)
        ea_mkdir([directory,'fiducials',filesep,'native']);
        ea_mkdir([directory,'fiducials',filesep,ea_getspace]);
        if exist([ea_space,'fiducials',filesep,fids(fid).name],'file')  % pair exists
            movefile([ea_space,'fiducials',filesep,fids(fid).name],...
                [directory,'fiducials',filesep,ea_getspace,filesep,fids(fid).name]);
            movefile([directory,'fiducials',filesep,fids(fid).name],...
                [directory,'fiducials',filesep,'native',filesep,fids(fid).name]);
        end
    end
    
end

