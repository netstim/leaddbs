function ea_genflippedjointnii(rightFile, leftFile)
% Helper function to flip right and left VTA/E-field

files = {rightFile; leftFile};

for f=1:2
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

    if ~isfile(flippedFile)
        ea_flip_lr_nonlinear(files{f}, flippedFile, 0);
    end
end
