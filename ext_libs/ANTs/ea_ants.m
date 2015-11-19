function ea_ants(fixedimage, movingimage, outputname)
% Wrapper for ANTs

if fileparts(movingimage)
    volumedir = [fileparts(movingimage), filesep];
else
    volumedir =['.', filesep];
end

ttries=ea_detct2anatattempts(volumedir); % how many coregistration attempts have been made in the past..

% first attempt..
fixparams{1} =[' 3' ...
             ' -t affine[ 0.1 ]' ...
             ' --regularization Gauss[3,0]' ...
             ' --number-of-iterations 10x10x5' ...
             ' --rigid-affine true'...
             ' --do-rigid true'];

% second attempt..
fixparams{2} = [' 3' ...
             ' -t affine[ 0.1 ]' ...
             ' --regularization Gauss[3,0]' ...
             ' --number-of-iterations 10x10x5' ...
             ' --initial-affine ', volumedir, 'ct2anat',num2str(ttries),'.txt'...
             ' --rigid-affine true'...
             ' --do-rigid true'];

if ttries>2
    ttries=2;
end


basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    ANTS = [basedir, 'ANTS.exe'];
    Warp = [basedir, 'WarpImageMultiTransform.exe'];
    
elseif isunix
    ANTS = [basedir, 'ANTS.', computer];
    Warp = [basedir, 'WarpImageMultiTransform.', computer];
end

ea_libs_helper
system([ANTS, fixparams{ttries}, ...
        ' -m PR[', fixedimage,',', movingimage,',1,2]' ...
        ' --output-naming ', outputname]);

outputbase = outputname(1:end-4);
system([Warp ...
        ' 3 ', movingimage, ' ', outputname, ...
        ' -R ', fixedimage, ...
        ' ', outputbase, 'Warp.nii' ...
        ' ', outputbase, 'Affine.txt']);

delete([outputbase, 'Warp.nii']);
delete([outputbase, 'InverseWarp.nii']);
movefile([outputbase, 'Affine.txt'], [volumedir, 'ct2anat.txt']);