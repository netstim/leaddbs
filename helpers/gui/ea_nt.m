function estring=ea_nt(options)

switch options.native
    case 1
        estring=['native',filesep];
    case 0
        estring=[ea_getspace,filesep];
end
