function ea_ants(fixedimage, movingimage, outputname)
% Wrapper for ANTs

fixparams = [' 3' ...
             ' --transformation-model SyN[0.3]' ...
             ' --regularization Gauss[3,0]' ...
             ' --number-of-iterations 10x10x5'];
     
if fileparts(movingimage)
    imagedir = [fileparts(movingimage), filesep];
else
    imagedir =['.', filesep];
end

if exist([imagedir, 'ct2anat.txt'],'file')
   fixparams = [fixparams, ' --initial-affine ', imagedir, 'ct2anat.txt']; 
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
system([ANTS, fixparams, ...
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
movefile([outputbase, 'Affine.txt'], [imagedir, 'ct2anat.txt']);