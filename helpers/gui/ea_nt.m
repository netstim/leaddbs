function estring = ea_nt(options)
% Return the space name + filesep (to construct folder path)

if isstruct(options)
    switch options.native
        case 1
            estring=['native',filesep];
        case 0
            estring=[ea_getspace,filesep];
    end
elseif isnumeric(options)
    switch options
        case 1
            estring=['native',filesep];
        case 0
            estring=[ea_getspace,filesep];
    end
elseif ischar(options)
    switch options
        case 'native'
            estring=['native',filesep];
        case 'mni'
            estring=[ea_getspace,filesep];
    end
end
