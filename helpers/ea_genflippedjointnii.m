function ea_genflippedjointnii(rightFile, leftFile)
% Helper function to flip right and left VTA/E-field

files = {rightFile; leftFile};

for f=1:2
    if ~isfile(files{f})
        ea_cprintf('CmdWinWarnings', 'Skipping flipping image: VTA doesn''t exist!\n');
        continue;
    end

    [fPath, fName, fExt] = fileparts(files{f});
    if ~isBIDSFileName(files{f})
        flippedFile = fullfile(fPath, ['fl_', fName, fExt]);
    else
        switch f
            case 1
                flippedFile = setBIDSEntity(files{f}, 'hemi', 'L', 'hemidesc', 'FlippedFromRight');
            case 2
                flippedFile = setBIDSEntity(files{f}, 'hemi', 'R',  'hemidesc', 'FlippedFromLeft');
        end
    end

    if contains(flippedFile, 'efield')
        interp = 1;
    else
        interp = 0;
    end

    if ~isfile(flippedFile)
        ea_flip_lr_nonlinear(files{f}, flippedFile, interp);
    end
end
