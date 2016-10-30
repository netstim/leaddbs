function ea_flirt(varargin)
% Wrapper for FSL linear registration

fixedimage=varargin{1};
movingimage=varargin{2};
outputimage=varargin{3};

if nargin>3
    writematout=varargin{4};
else
    writematout=1;
end
if nargin>4
    otherfiles=varargin{5};
end

if fileparts(movingimage)
    volumedir = [fileparts(movingimage), filesep];
else
    volumedir =['.', filesep];
end

[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

fixedimage = ea_path_helper(ea_niigz(fixedimage));
movingimage = ea_path_helper(ea_niigz(movingimage));
outputimage = ea_path_helper(ea_niigz(outputimage));

basedir = [fileparts(mfilename('fullpath')), filesep];
setenv('FSLOUTPUTTYPE','NIFTI')
if ispc
    FLIRT = [basedir, 'flirt.exe'];
else
    FLIRT = [basedir, 'flirt.', computer('arch')];
end


% determine how many runs have been performed before..
runs=0;
for r=1:5;
    if exist([volumedir, 'fslaffine',num2str(r),'.mat'],'file')
        runs=r;
    else
        break
    end
end

if runs==0 % mattes MI affine + rigid

    affinestage = [' -cost mutualinfo' ...
    ' -searchcost mutualinfo' ...
    ' -interp sinc' ...
    ' -omat fslaffine',num2str(runs),'.mat' ...
    ' -verbose 1'];

elseif runs>0

    affinestage = [' -init fslaffine',num2str(runs),'.mat' ...
    ' -cost mutualinfo' ...
    ' -searchcost mutualinfo' ...
    ' -interp sinc' ...
    ' -omat ct2anat',num2str(runs),'.mat' ...
    ' -verbose 1'];

end

ea_libs_helper
cmd = [FLIRT, ...
    ' -ref ',fixedimage, ...
    ' -in ',movingimage, ...
    ' -out ',outputimage ...
    affinestage];
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

if exist('otherfiles','var')
    ea_error('This procedure is not yet supported using FSL. Please choose a different coregistration method.');
    if ~isempty(otherfiles)
        for ofi=1:length(otherfiles)
        [options.root,options.patientname]=fileparts(fileparts(otherfiles{ofi}));
        options.root=[options.root,filesep];
        options.prefs=ea_prefs(options.patientname);
        ea_ants_applytransforms(options,otherfiles(ofi),otherfiles(ofi),0,fixedimage,[outputbase, '0GenericAffine.mat']);
        end
    end
end


