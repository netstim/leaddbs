function ea_fnirt(varargin)
% Wrapper for FSL FNIRT nonlinear registration

fixedimage=varargin{1};
movingimage=varargin{2};
outputimage=varargin{3};


[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

if ischar(fixedimage)
    fixedimage={fixedimage};
elseif ~iscell(fixedimage)
	ea_error('Please supply variable fixedimage as either char or cellstring');
end



if nargin>3
    weights=varargin{4};
    metrics=varargin{5};
%     options=varargin{6};
else
    weights=ones(length(fixedimage),1);
    metrics=repmat({'MI'},length(fixedimage),1);
end

if ischar(movingimage)
    movingimage={movingimage};
elseif ~iscell(movingimage)
    ea_error('Please supply variable fixedimage as either char or cellstring');
end

directory=fileparts(movingimage{1});
directory=ea_path_helper([directory,filesep]);



for fi=1:length(fixedimage)
    fixedimage{fi} = ea_path_helper(ea_niigz(fixedimage{fi}));
end
for fi=1:length(movingimage)
    movingimage{fi} = ea_path_helper(ea_niigz(movingimage{fi}));
end

if length(fixedimage)~=length(movingimage)
    ea_error('Please supply pairs of moving and fixed images (can be repetitive).');
end

outputimage = ea_path_helper(ea_niigz(outputimage));

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    FLIRT = [basedir, 'flirt.exe'];
    BET = [basedir, 'bet2.exe'];
    FNIRT = [basedir, 'fnirt.exe'];
    INVWARP = [basedir, 'invwarp.exe'];
else
    FLIRT = [basedir, 'flirt.', computer('arch')];    
    BET = [basedir, 'bet2.', computer('arch')];    
    FNIRT = [basedir, 'fnirt.', computer('arch')];
    INVWARP = [basedir, 'invwarp.', computer('arch')];
end


setenv('FSLOUTPUTTYPE','NIFTI')
flirtstage = [' -ref ',fixedimage{1},' -in ',movingimage{1},' -omat ',directory,'aff_anat2mni.mat -verbose 1 -cost mutualinfo -searchcost mutualinfo'];
betstage = [' ',movingimage{1},' ',ea_prependb(movingimage{1}),' -m'];
[bfile,mask,weirdmask]=ea_prependb(movingimage{1});
prefs=ea_prefs('');
%fnirtstage = [' --ref=',fixedimage{1},' --in=',movingimage{1},' --aff=',directory,'aff_anat2mni.mat --inmask=',ea_path_helper(mask),' --refmask=',ea_path_helper([ea_getearoot,'templates',filesep,'bmni_hires.nii.gz']),' --iout=',outputimage,' --cout=',directory,'warp_struct2mni.nii --verbose --subsamp=4,4,2,2,1,1 --warpres=8,8,8 --miter=5,5,5,5,5,10 --infwhm=8,6,5,4.5,3,2 --reffwhm=8,6,5,4,2,0 --lambda=300,150,100,50,40,30 --estint=1,1,1,1,1,0 --ssqlambda=1 --regmod=bending_energy --intmod=global_non_linear_with_bias --intorder=5 --biasres=50,50,50 --biaslambda=10000 --refderiv=0 --applyrefmask=1,1,1,1,1,1 --applyinmask=1'];
fnirtstage = [' --ref=',fixedimage{1},' --in=',movingimage{1},' --aff=',directory,'aff_anat2mni.mat --inmask=',ea_path_helper(mask),' --refmask=',ea_path_helper([ea_getearoot,'templates',filesep,'bmni_hires.nii.gz']),' --iout=',outputimage,' --cout=',directory,prefs.gprenii,'Composite.nii --verbose --subsamp=4,2 --warpres=8,8 --miter=5,10 --infwhm=8,2 --reffwhm=4,0 --lambda=300,30 --estint=1,0 --ssqlambda=1 --regmod=bending_energy --intmod=global_non_linear_with_bias --intorder=5 --biasres=50,50,50 --biaslambda=10000 --refderiv=0 --applyrefmask=1,1 --applyinmask=1,1'];

invwarpstage = [' -w ',directory,'glanatComposite.nii -o ',directory,prefs.gprenii,'InverseComposite.nii -r ',movingimage{1}];

ea_libs_helper

lcmd = [FLIRT, flirtstage];
bcmd = [BET, betstage];
ncmd = [FNIRT, fnirtstage];
icmd = [INVWARP, invwarpstage];

if ~ispc
    system(['bash -c "', lcmd, '"']);    
    system(['bash -c "', bcmd, '"']);  
    ea_cleanupmask(movingimage{1})
    system(['bash -c "', ncmd, '"']);
    system(['bash -c "', icmd, '"']);
else
    system(lcmd);    
    system(bcmd);
    ea_cleanupmask(movingimage{1})
    system(ncmd);
    system(icmd);
end

function [fname,maskfn,weirdmaskfn]=ea_prependb(fname)
[pth,fn,ext]=fileparts(fname);
fname=fullfile(pth,['b',fn,ext]);
maskfn=fullfile(pth,['b',fn,'_mask',ext]);
weirdmaskfn=fullfile(pth,['b',fn,'.nii_mask',ext]);

function ea_cleanupmask(fname)   
[bfile,mask,weirdmask]=ea_prependb(fname);
movefile(weirdmask,mask)
delete(bfile)