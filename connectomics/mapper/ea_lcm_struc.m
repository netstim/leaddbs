function ea_lcm_struc(options)


if strcmp(options.lcm.struc.connectome,'No structural connectome found.')
    return
end
disp('Running structural connectivity...');

if isempty(options.uipatdirs)
    base=ea_getconnectomebase();
else
    base=[options.uipatdirs,'connectomes',filesep,'data',filesep];
end
cs_dmri_conseed(base,options.lcm.struc.connectome,...
    options.lcm.seeds,...
    ea_lcm_resolvecmd(options.lcm.cmd),...
    '0',...
    options.lcm.odir,...
    options.lcm.omask,...
    ea_resolve_espace(options.lcm.struc.espace));
disp('Done.');

function fi=ea_resolve_espace(sp)

base=ea_getconnectomebase;
switch sp
    case 1
        fi=['222.nii'];
    case 2
        fi=['111.nii'];
    case 3
        fi=['555.nii'];
end