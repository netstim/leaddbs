function [cmdout, tempf] = ea_dcm2nii(inputimage, outputimage)
% Wrapper for dcm2nii, just used for reorientation and cropping currently

% Save the result to a new file, when the second parameter is specified
if nargin == 2
    copyfile(inputimage,outputimage)
    inputimage = outputimage;
end

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    dcm2nii = [basedir, 'dcm2nii.exe'];
else
    dcm2nii = [basedir, 'dcm2nii.', computer('arch')];
end

cmd=[dcm2nii, ' -g n -x y ', inputimage];

fprintf('\nReorient and crop image...\n\n');
if ~ispc
    [~,cmdout] = system(['bash -c "', cmd, '"']);
else
    [~,cmdout] = system(cmd);
end

disp(cmdout);

% Check the output files of dcm2nii
savedf = regexp(cmdout, '(?<=Saving )\S+','match');

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
