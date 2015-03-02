function varargout=ea_normalize_apply_normalization(options)
% This is a function that simply applies normalization parameters (y_ea_normparams.nii) to the
% unnormalized files.
%
% Since a high resolution is needed for accurate DBS localizations, this
% function applies DARTEL to an output resolution of 0.5 mm isotropic. This
% makes the procedure quite slow.

% The function uses some code snippets written by Ged Ridgway.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='(Re-)apply (priorly) estimated normalization.';
    varargout{2}={'SPM8','SPM12'};
    return
end


% Apply estimated deformation to (coregistered) post-op data.
switch spm('ver')
    
    case 'SPM8'
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_normparams.nii']};
        matlabbatch{1}.spm.util.defs.ofname = '';
        
        postops={options.prefs.tranii_unnormalized,options.prefs.cornii_unnormalized,options.prefs.sagnii_unnormalized,options.prefs.prenii_unnormalized,options.prefs.ctnii_coregistered};
        cnt=1;
        for postop=1:length(postops)
            if exist([options.root,options.patientname,filesep,postops{postop}],'file')
                matlabbatch{1}.spm.util.defs.fnames{cnt}=[options.root,options.patientname,filesep,postops{postop},',1'];
                cnt=cnt+1;
            end
        end
        
        matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.patientname,filesep]};
        matlabbatch{1}.spm.util.defs.interp = 1;
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear matlabbatch jobs;
        
        % rename files:
        try copyfile([options.root,options.patientname,filesep,'w',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gprenii]); end
        try movefile([options.root,options.patientname,filesep,'w',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,options.prefs.prenii]); end
        try copyfile([options.root,options.patientname,filesep,'w',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gtranii]); end
        try movefile([options.root,options.patientname,filesep,'w',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.tranii]); end
        try copyfile([options.root,options.patientname,filesep,'w',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gcornii]); end
        try movefile([options.root,options.patientname,filesep,'w',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.cornii]); end
        try copyfile([options.root,options.patientname,filesep,'w',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gsagnii]); end
        try movefile([options.root,options.patientname,filesep,'w',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.sagnii]); end
        try copyfile([options.root,options.patientname,filesep,'w',options.prefs.ctnii_coregistered],[options.root,options.patientname,filesep,options.prefs.gctnii]); end
        try movefile([options.root,options.patientname,filesep,'w',options.prefs.ctnii_coregistered],[options.root,options.patientname,filesep,options.prefs.ctnii]); end
    case 'SPM12'
        % export lfiles (fine resolution, small bounding box.
        postops={options.prefs.tranii_unnormalized,options.prefs.cornii_unnormalized,options.prefs.sagnii_unnormalized,options.prefs.prenii_unnormalized,options.prefs.ctnii_coregistered};
        
        for postop=1:length(postops)
            if exist([options.root,options.patientname,filesep,postops{postop}],'file')
                nii=ea_load_untouch_nii([options.root,options.patientname,filesep,postops{postop}]);
                gaussdim=abs(nii.hdr.dime.pixdim(2:4));
                resize_img([options.root,options.patientname,filesep,postops{postop}],gaussdim./2,nan(2,3),0);
                %gaussdim=abs(gaussdim(1:3)).*2;
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.push.fnames{1}=[options.root,options.patientname,filesep,'r',postops{postop},''];
                matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
                matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {[options.root,options.patientname,filesep]};
                matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.bb = [-55 45 9.5; 55 -65 -25];
                matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.vox = [0.22 0.22 0.22];
                matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
                matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = gaussdim;
                jobs{1}=matlabbatch;
                cfg_util('run',jobs);
                clear matlabbatch jobs;
            end
        end
        
        
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,options.prefs.prenii]); end
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.tranii]); end
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.cornii]); end
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.sagnii]); end
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.ctnii_coregistered],[options.root,options.patientname,filesep,options.prefs.ctnii]); end
        
        % export glfiles (a bit more coarse resolution, full brain bounding box.
        
        for postop=1:length(postops)
            if exist([options.root,options.patientname,filesep,postops{postop}],'file')
                nii=ea_load_untouch_nii([options.root,options.patientname,filesep,postops{postop}]);
                gaussdim=abs(nii.hdr.dime.pixdim(2:4));
               %gaussdim=abs(gaussdim(1:3)).*2;
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.push.fnames{1}=[options.root,options.patientname,filesep,'r',postops{postop},''];
                matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
                matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {[options.root,options.patientname,filesep]};
                matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.bb = [-78 -112  -50
                    78   76   85];
                matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.vox = [0.5 0.5 0.5];
                matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
                matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = gaussdim;
                jobs{1}=matlabbatch;
                cfg_util('run',jobs);
                clear matlabbatch jobs;
            end
        end
        
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gprenii]); end
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gtranii]); end
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gcornii]); end
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.gsagnii]); end
        try movefile([options.root,options.patientname,filesep,'swr',options.prefs.ctnii_coregistered],[options.root,options.patientname,filesep,options.prefs.gctnii]); end
                try delete([options.root,options.patientname,filesep,'r',options.prefs.prenii_unnormalized]); end
        try delete([options.root,options.patientname,filesep,'r',options.prefs.tranii_unnormalized]); end
        try delete([options.root,options.patientname,filesep,'r',options.prefs.cornii_unnormalized]); end
        try delete([options.root,options.patientname,filesep,'r',options.prefs.sagnii_unnormalized]); end
        try delete([options.root,options.patientname,filesep,'r',options.prefs.ctnii_coregistered]); end
        
end


function [flirtmat spmvoxmat fslvoxmat] = worldmat2flirtmat(worldmat, src, trg)
%worldmat2flirtmat: convert NIfTI world (mm) coordinates matrix to flirt
%
% Example:
%  [flirtmat spmvoxmat fslvoxmat] = worldmat2flirtdmat(worldmat, src, trg);
%
% See also: flirtmat2worldmat, flirtmat_write

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

if ischar(src)
    src = nifti(src);
end
if ischar(trg)
    trg = nifti(trg);
end

spmvoxmat = inv(src.mat) * worldmat * trg.mat;
addone = eye(4); addone(:, 4) = 1;
fslvoxmat = inv(addone) * spmvoxmat * addone;
trgscl = nifti2scl(trg);
srcscl = nifti2scl(src);
flirtmat = inv( srcscl * fslvoxmat * inv(trgscl) );



function flirtmat_write(fname, mat)
%flirtmat_write: save a 4-by-4 matrix to a file as handled by flirt -init
% Example:
%  flirtmat_write(fname, affinemat)
% See also: flirtmat_read, flirtmat2worldmat, worldmat2flirtmat

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

fid = fopen(fname, 'w');
% Note that transpose is needed to go from MATLAB's column-major matrix to
% writing the file in a row-major way (see also flirtmat_read)
str = sprintf('%f  %f  %f  %f\n', mat');
fprintf(fid, str(1:end-1)); % drop final newline
fclose(fid);




function [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg)
%flirtmat2worldmat: convert saved flirt matrix to NIfTI world coords matrix
% flirt matrix is from text file specified in "flirt -omat mat.txt" command
% world matrix maps from NIfTI world coordinates in target to source. Note:
% mat.txt contains a mapping from source to target in FSL's *scaled* coords
% which are not NIfTI world coordinates, and note src-trg directionality!
% worldmat from this script reproduces "img2imgcoord -mm ...".
%
% The script can also return a matrix to map from target to source voxels
% in MATLAB/SPM's one-based convention, or in FSL's zero-based convention
%
% Example:
%  [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
%
% See also: worldmat2flirtmat, flirtmat_read, flirtmat_write

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

if ischar(src)
    src = nifti(src);
end
if ischar(trg)
    trg = nifti(trg);
end

if ischar(flirtmat)
    flirtmat = flirtmat_read(flirtmat);
end

% src = inv(flirtmat) * trg
% srcvox = src.mat \ inv(flirtmat) * trg.mat * trgvox
% BUT, flirt doesn't use src.mat, only absolute values of the
% scaling elements from it,
% AND, if images are not radiological, the x-axis is flipped, see:
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0810&L=FSL&P=185638
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0903&L=FSL&P=R93775
trgscl = nifti2scl(trg);
srcscl = nifti2scl(src);
fslvoxmat = inv(srcscl) * inv(flirtmat) * trgscl;

% AND, Flirt's voxels are zero-based, while SPM's are one-based...
addone = eye(4); addone(:, 4) = 1;
spmvoxmat = addone * fslvoxmat * inv(addone);

worldmat = src.mat * spmvoxmat * inv(trg.mat);


%%
function scl = nifti2scl(N)
% not sure if this is always correct with rotations in mat, but seems okay!
scl = diag([sqrt(sum(N.mat(1:3,1:3).^2)) 1]);
if det(N.mat) > 0
    % neurological, x-axis is flipped, such that [3 2 1 0] and [0 1 2 3]
    % have the same *scaled* coordinates:
    xflip = diag([-1 1 1 1]);
    xflip(1, 4) = N.dat.dim(1)-1; % reflect about centre
    scl = scl * xflip;
end

