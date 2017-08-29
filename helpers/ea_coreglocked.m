function locked=ea_coreglocked(options,moving)
locked=0;
if options.overwriteapproved
    return
end
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'ea_coreg_approved.mat'],'file')
    return
else
    approved=load([directory,'ea_coreg_approved.mat']);
    if isfield(approved,striponeex(moving))
        if approved.(striponeex(moving))
            locked=approved.(striponeex(moving));
        end
    end 
end

function str=striponeex(str)

str=strrep(str,',1','');

str=strrep(str,'.nii','');

str=strrep(str,'.gz','');
[~,str]=fileparts(str);
