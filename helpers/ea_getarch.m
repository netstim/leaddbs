function arch = ea_getarch
% Return computer arch, hanlding Mac arch properly for MATLAB ver < 2023b
% on Apple Silicon.

if ~ismac
    arch = computer('arch');
else
    if isMatlabVer('<', [23,2])
        [~, cmdout] = system('uname -v');
        if contains(cmdout, 'ARM64')
            arch = 'maca64';
        else
            arch = 'maci64';
        end
    else
        arch = computer('arch');
    end
end
