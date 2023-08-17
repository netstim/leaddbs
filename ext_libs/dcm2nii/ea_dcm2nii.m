function [cmdout, tempf] = ea_dcm2nii(inputimage, outputimage)
% Wrapper for dcm2nii, only used for reorientation and cropping
%
% dcm2ii is not stable enough in some cases (will write out corrupted file
% with 1KB size). Use 'ea_rocrop' instead to do reorientation and cropping.

% Save to a new file, when 'outputimage' is specified.
if nargin == 2
    copyfile(inputimage, outputimage)
    inputimage = outputimage;
end

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

dcm2nii = ea_getExec([basedir, 'dcm2nii'], escapePath = 1);


cmd=[dcm2nii, ' -g n -x y ', ea_path_helper(inputimage)];

fprintf('\nReorient and crop image...\n\n');
[~, cmdout] = ea_runcmd(cmd);

disp(cmdout);

if contains(cmdout,'Crop 3D unsupported datatype.') ||...
   contains(cmdout,'Serious error: format not supported by Crop3D.') ||...
   contains(cmdout,'niftiutil UnSwapImg error: datatype not supported.')

    fprintf(['Cropping failed due to improper datatype...\nTrying to ' ...
             'change the datatype...\n']);

    % Change datatype
    nii = load_untouch_nii(inputimage);
    switch nii.hdr.dime.datatype
        case 512
            nii.hdr.dime.datatype = 4;
            fprintf('datatype changed from UINT16 to INT16...\n');
        case 768
            nii.hdr.dime.datatype = 8;
            fprintf('datatype changed from UINT32 to INT32...\n')
        % otherwise
        %     nii.hdr.dime.datatype = 16;
        %     fprintf('datatype changed to SINGLE/FLOAT32...\n')
    end
    save_untouch_nii(nii, inputimage);

    % Try again now
    fprintf('\nRetry reorienting and cropping...\n\n');
    [~, cmdout] = ea_runcmd(cmd);

    disp(cmdout);
end

% Check the output files of dcm2nii
cmdout = strsplit(cmdout, '\x{0A}', 'DelimiterType', 'RegularExpression')';
savedf = regexp(cmdout, '(?<=Saving )(.*?)$', 'match', 'once');
savedf(cellfun(@isempty, savedf)) = [];

if ~isempty(savedf)
    movefile(savedf{end}, inputimage);
    if numel(savedf) == 2
        delete(savedf{1});
    end
    tempf = savedf{end};
    disp('Reorientation and/or cropping applied.');
else
    tempf = inputimage;
    disp('No need to reorient or crop!');
end
