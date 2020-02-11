function ea_lcm_struc(options)

if strcmp(options.lcm.struc.connectome,'No structural connectome found.')
    return
end
disp('Running structural connectivity...');

if strcmp(options.lcm.struc.connectome, 'Patient''s fiber tracts')
    base='Patient''s fiber tracts';
    options.lcm.struc.connectome=options.prefs.FTR_normalized;
else
    base=ea_getconnectomebase();
end

cs_dmri_conseed(base,options.lcm.struc.connectome,...
    options.lcm.seeds',...
    ea_lcm_resolvecmd(options.lcm.cmd),...
    options.lcm.odir,...
    ea_resolve_espace(options.lcm.struc.espace),options);
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
