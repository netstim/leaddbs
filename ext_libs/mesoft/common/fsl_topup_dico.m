function fsl_topup_dico(reversedPE_b0File, targetFiles, varargin)

if nargin >2
    tempDirName = varargin{1};
else
    tempDirName = 'fsl_tempDir';
end

%if the function is called via spm: remove these numbers for 4d data at the
% end
if ~strcmp(reversedPE_b0File(end-2:end), 'nii')
    reversedPE_b0File = reversedPE_b0File(1:end-2);
    targetFiles = cellfun( @(x) x(1:end-2), targetFiles, 'uniformoutput', 0');
end


%% set FSL PATH if not set yet
if isempty(strfind(getenv('LD_LIBRARY_PATH'), '/fsl/'))
    setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') ':/usr/lib/fsl/5.0']);
    setenv('PATH', [getenv('PATH') ':/usr/lib/fsl/5.0']);
    setenv('FSLOUTPUTTYPE', 'NIFTI');
end

%% find filenames and paths
if exist(tempDirName, 'dir')
    system(['rm -r ' tempDirName]);
end
mkdir(tempDirName);

for k = 1:length(targetFiles)
    tFile = targetFiles{k};
    system(['cp ' tFile ' ' tempDirName]);
    system(['cp ' [tFile(1:end-4) '.mat'] ' ' tempDirName]);  % must also copy the mat files for diffusion
end

fileList  = targetFiles;
PAb0Name = reversedPE_b0File;
APb0Name = targetFiles{1};

% define tempoorary files
b02b0       = fullfile(tempDirName, 'b02b0.cnf');
acqparams   = fullfile(tempDirName, 'acqparams.txt');
AP_PA_b0    = fullfile(tempDirName, 'AP_PA_b0.nii');
AP_PA_b0Extended= fullfile(tempDirName, 'AP_PA_b0Extended.nii');
myTempFile  = fullfile(tempDirName, 'myTempFile.nii');
tempSlice1  = fullfile(tempDirName, 'tempSlice1.nii');
tempSlice2  = fullfile(tempDirName, 'tempSlice2.nii');
topup_AP_PA_b0 = fullfile(tempDirName, 'topup_AP_PA_b0');



%%  merge AP and PA images into timeseries of images
system(['fslmerge -t  ' AP_PA_b0 ' ' APb0Name ' ' PAb0Name ]);

% load image in order to get image parameterd
tImg = nifti(AP_PA_b0);
tSize = size(tImg.dat);



%% extend image volume
% extract first and last slice, and extend image. This must be done
% because of some circular effects of topup??
system(['fslroi ' AP_PA_b0 ' ' tempSlice1 ' 0 -1 0 -1 0 1 0 -1']);
system(['fslroi ' AP_PA_b0 ' ' tempSlice2 ' 0 -1 0 -1 ' num2str(tSize(3)-1) ' 1 0 -1']);
system(['fslmerge -z ' AP_PA_b0Extended ' ' tempSlice1 ' ' AP_PA_b0 ' ' tempSlice2 ]);
system(['fslmerge -z ' AP_PA_b0Extended ' ' tempSlice1 ' ' tempSlice1 ' ' tempSlice1 ' ' AP_PA_b0 ' ' tempSlice2 ]);
%system(['fslmerge -z ' AP_PA_b0Extended ' ' tempSlice1 ' ' tempSlice1 ' ' tempSlice1 ' ' tempSlice1 ' ' tempSlice1 ' ' tempSlice1 ' ' AP_PA_b0 ' ' tempSlice2 ]);

% add additional slice if nSlices is odd, otherwis subsampling will not
% match.
if mod(tSize(3),2)
    system(['fslmerge -z ' AP_PA_b0Extended ' ' AP_PA_b0Extended ' ' tempSlice2]);
end


%%  configure and run the topup 
% write the directions into a file
dlmwrite(acqparams, [0  1  0 0.0665;     0 -1 0 0.0665], 'delimiter', ' ');


fid = fopen(b02b0, 'w+');
params = '--warpres=20,16,14,12,10,6,4,4,4 ';fprintf(fid, '%s\n', params);
% Subsampling level (a value of 2 indicates that a 2x2x2 neighbourhood is collapsed to 1 voxel)
params = '--subsamp=2,2,2,2,2,1,1,1,1 '; fprintf(fid, '%s\n', params);
% FWHM of gaussian smoothing
params = '--fwhm=8,6,4,3,3,2,1,0,0 ';fprintf(fid, '%s\n', params);
% Maximum number of iterations
params = '--miter=5,5,5,5,5,10,10,20,20 ';fprintf(fid, '%s\n', params);
% Relative weight of regularisation
params = '--lambda=0.005,0.001,0.0001,0.000015,0.000005,0.0000005,0.00000005,0.0000000005,0.00000000001 ';fprintf(fid, '%s\n', params);
% If set to 1 lambda is multiplied by the current average squared difference
params = '--ssqlambda=1 ';fprintf(fid, '%s\n', params);
% Regularisation model
params = '--regmod=bending_energy ';fprintf(fid, '%s\n', params);
% If set to 1 movements are estimated along with the field
params = '--estmov=0,0,0,0,0,0,0,0,0';fprintf(fid, '%s\n', params);
% 0=Levenberg-Marquardt, 1=Scaled Conjugate Gradient
params = '--minmet=0,0,0,0,0,1,1,1,1 ';fprintf(fid, '%s\n', params);
% Quadratic or cubic splines
params = '--splineorder=3 ';fprintf(fid, '%s\n', params);
% Precision for calculation and storage of Hessian
params = '--numprec=double ';fprintf(fid, '%s\n', params);
% Linear or spline interpolation
params = '--interp=spline ';fprintf(fid, '%s\n', params);
% If set to 1 the images are individually scaled to a common mean intensity 
params = '--scale=1 ';fprintf(fid, '%s\n', params);
params = '--verbose ';fprintf(fid, '%s\n', params);
fclose(fid);

% ooold system('topup --imain=' AP_PA_b0Extended ' --datain=' acqparams ' --config=' b02b0.cnf ' --out=' topup_AP_PA_b0 ' --fout=my_field --iout=my_unwarped_images');
%system(['topup --imain=' AP_PA_b0Extended ' --datain=' acqparams ' --config=' b02b0 ' --out=' topup_AP_PA_b0]);

system(['topup --imain=' AP_PA_b0Extended ' --datain=' acqparams ' --config=' b02b0 ' --out=' topup_AP_PA_b0 ' --fout=my_field --iout=my_unwarped_images']);


%% apply the warp field to all images.

for k=1:length(fileList);
    tFile = fileList{k};
    disp(['Applying topup to ' tFile])
    
    system(['fslroi ' tFile  ' ' tempSlice1 ' 0 -1 0 -1 0 1 0 -1']);
    system(['fslroi ' tFile  ' ' tempSlice2 ' 0 -1 0 -1 ' num2str(tSize(3)-1) ' 1 0 -1']);
    system(['fslmerge -z ' myTempFile ' ' tempSlice1 ' ' tFile ' ' tempSlice2 ]);
    
    % add additional slice if nSlices is odd, otherwis subsampling will not
    % match.
    if mod(tSize(3),2)
        system(['fslmerge -z '  myTempFile ' ' myTempFile ' ' tempSlice2]);
    end
    
    [tPath,tFileName,tExt] = fileparts(tFile);
    tFileOut = fullfile(tempDirName, [tFileName tExt]);
    
    system(['applytopup --imain=' myTempFile  ' --out=' tFileOut '  --topup=' topup_AP_PA_b0 ' --datain=' acqparams ' --inindex=1 --verbose --method=jac']);
    system(['fslroi ' tFileOut   ' ' tFileOut ' 0 -1 0 -1  1 ' num2str(tSize(3)) '  0 -1']);
    
end


%% cleanup
% system('rm b02b0.cnf');
% system('rm acqparams.txt');
% 
% system(['rm tempSlice1.nii' ]);
% system(['rm tempSlice2.nii' ]);
% system(['rm AP_PA_b0.nii' ]);
% system(['rm myTempFile.nii' ]);
% 
% system(['rm ' b02b0]);
% system(['rm ' acqparams]);
% system(['rm ' tempSlice1 ]);
% system(['rm ' tempSlice2 ]);
% system(['rm ' AP_PA_b0 ]);
% system(['rm ' AP_PA_b0Extended '*'])
% system(['rm ' myTempFile ]);
% system(['rm ' topup_AP_PA_b0  '*']);


%system(['rm my_field.nii' ])
%system(['rm my_unwarped_images.nii' ])
%system(['rm topup_AP_PA_b0_fieldcoef.nii' ])
%system(['rm topup_AP_PA_b0_movpar.txt' ])


