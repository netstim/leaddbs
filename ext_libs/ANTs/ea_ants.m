function ea_ants(fixedimage, movingimage, outputname)
% Wrapper for ANTs nonlinear registration

if fileparts(movingimage)
    volumedir = [fileparts(movingimage), filesep];
else
    volumedir =['.', filesep];
end

ttries=ea_detct2anatattempts(volumedir); % how many coregistration attempts have been made in the past..

% first attempt..
fixparams{1} =[' -d 3 -t a -v 1 -c [10000x0x0,1.e-8,20] -s 4x2x1vox'];

fixparams{1} =[' -d 3 -v 1 -t translation[ 0.1 ]  -c [10,1.e-8,20]   -s 4x2x1vox  -f 6x4x2 -l 1  -t rigid[ 0.1 ] -c [10,1.e-8,20]  -s 4x2x1vox   -f 3x2x1 -l 1  -t affine[ 0.1 ]   -c [10,1.e-8,20] -s 4x2x1vox']


% second attempt..
%fixparams{2} = [' -d 3 -t a -v 1 -c 1e-6 -s 1e-6'];

if ttries>1
    ttries=1;
end


basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    ANTS = [basedir, 'antsRegistration.exe'];

elseif isunix
    ANTS = [basedir, 'antsRegistration.', computer];
end

ea_libs_helper
system([ANTS,' -m ',movingimage,'.gz',' -f ',fixedimage,'.gz', ' -o ',outputname,fixparams{ttries}]);


    keyboard
movefile([outputbase, 'Affine.txt'], [volumedir, 'ct2anat', num2str(ttries+1), '.txt']);
