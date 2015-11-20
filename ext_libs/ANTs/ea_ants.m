function ea_ants(fixedimage, movingimage, outputname)
% Wrapper for ANTs nonlinear registration

if fileparts(movingimage)
    volumedir = [fileparts(movingimage), filesep];
else
    volumedir =['.', filesep];
end

ttries=ea_detct2anatattempts(volumedir); % how many coregistration attempts have been made in the past..

% first attempt..
fixparams{1} =[' -d 3 -t a'];

% second attempt..
fixparams{2} = [' -d 3 -t a'];

if ttries>2
    ttries=2;
end


basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    ANTS = [basedir, 'antsRegistration.exe'];

elseif isunix
    ANTS = [basedir, 'antsRegistration.', computer];
end

ea_libs_helper
system([ANTS,' -m ',movingimage,' -f ',fixedimage, ' -o ',outputname,fixparams{ttries}]);


    keyboard
movefile([outputbase, 'Affine.txt'], [volumedir, 'ct2anat', num2str(ttries+1), '.txt']);
