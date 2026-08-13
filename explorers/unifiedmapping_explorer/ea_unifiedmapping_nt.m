function estring = ea_unifiedmapping_nt(calcspace)
% Return space folder using the unified-mapping convention:
%   0 = native
%   1 = template/MNI

switch calcspace
    case 0
        estring = ['native', filesep];

    case 1
        estring = [ea_getspace, filesep];

    otherwise
        error('ea_unifiedmapping_nt:InvalidCalcspace', ...
            'calcspace must be 0 (native) or 1 (template/MNI).');
end
end
