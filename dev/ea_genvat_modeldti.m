function varargout=ea_genvat_modeldti(varargin)
% This function generates a volume of activated tissue around for each
% electrode.
% Usage: VAT=ea_genvat(coords_mm,stimparams,options).
% ? stimparams is a struct variable with fields U (8*1 with voltage
% entries) and Im (8*1 with Impedance measurements).
%
% This function only touches the .VAT entry of stimparams struct of the
% given side.

if nargin==4
    coords=varargin{1};
    stimparams=varargin{2};
    side=varargin{3};
    options=varargin{4};
    elseif nargin==1

    if ischar(varargin{1}) % return name of method.
        varargout{1}='Model DTI';
        return
    end
end

%% build DTD (if needed)
if ~exist([options.root,options.patientname,filesep,options.prefs.DTD],'file')
    ea_prepare_dti(options);
end

if ~exist([options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_b0.mat'],'file')
disp('Building Mask in native DTI-space at electrode positions.');
%% generate mask

MNI=spm_vol([fileparts(which('spm')),filesep,'canonical',filesep,'single_subj_T1.nii']);
MNIX=spm_read_vols(MNI);
MNIX(:)=0;
% get indices of image
[xx,yy,zz]=ind2sub(size(MNIX),find(MNIX==0));
XYZ=[xx,yy,zz];
% transpose coords to voxel-space
coords_vox=[MNI.mat\[[coords{1};coords{2}],ones(size(coords{1},1)*2,1)]']';
coords_vox=coords_vox(:,1:3);
% get voxel ids around electrodes
idx=cell2mat(rangesearch(XYZ,coords_vox,10)');
MNIX(sub2ind(size(MNIX),XYZ(idx,1),XYZ(idx,2),XYZ(idx,3)))=1;

% save mask
mkdir([options.root,options.patientname,filesep,'vat_resources']);
MNI.fname=[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_mni.nii'];
spm_write_vol(MNI,MNIX);

% warp to subject space using inv parameters
switch spm('ver')
    case 'SPM8'
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']};
        matlabbatch{1}.spm.util.defs.ofname = '';
        matlabbatch{1}.spm.util.defs.fnames = {[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_mni.nii']};
        matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.patientname,filesep,'vat_resources']};
        matlabbatch{1}.spm.util.defs.interp = 1;
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
        clear jobs matlabbatch
        movefile([options.root,options.patientname,filesep,'vat_resources',filesep,'wtrackmask_mni.nii'],[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_pre.nii']);
    case 'SPM12'
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_normparams.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames = {[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_mni.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {[options.root,options.patientname,filesep,'vat_resources']};
        matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[options.root,options.patientname,filesep,options.prefs.tranii]};
        matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [4 4 4];
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
        clear jobs matlabbatch

        movefile([options.root,options.patientname,filesep,'vat_resources',filesep,'swtrackmask_mni.nii'],[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_pre.nii']);
end

copyfile([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'tmp.nii']);

% coreg to b0.
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[options.root,options.patientname,filesep,options.prefs.b0,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[options.root,options.patientname,filesep,'tmp.nii,1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_pre.nii']};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
jobs{1}=matlabbatch;
spm_jobman('run',jobs);
clear matlabbatch jobs

movefile([options.root,options.patientname,filesep,'vat_resources',filesep,'rtrackmask_pre.nii'],[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_b0.nii']);

% convert to maskstruct
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.srcimgs = {[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_b0.nii']};
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.thresh = [0.1 Inf];
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.nfval = false;
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.outchoice.outmat.outdir = {[options.root,options.patientname,filesep,'vat_resources']};
matlabbatch{1}.impexp_NiftiMrStruct.nifti2roistruct.outchoice.outmat.fname = 'trackmask_b0.mat';
jobs{1}=matlabbatch;
spm_jobman('run',jobs);
clear matlabbatch jobs

% clean up..
delete([options.root,options.patientname,filesep,'tmp.nii']);
delete([options.root,options.patientname,filesep,'rtmp.nii']);
delete([options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_mprage.nii']);
delete([options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_pre.nii']);
delete([options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_b0.nii']);
disp('Done building mask.');
end

%% track DTI
if ~exist([options.root,options.patientname,filesep,'vat_resources',filesep,'leadfield_FTR.mat'],'file')
    disp('Fibertracking around electrodes to generate leadfield, this will take a while...');
    matlabbatch{1}.dtijobs.tracking.GTtrack.fname.filenameHARDI = {[options.root,options.patientname,filesep,options.prefs.HARDI]};
    matlabbatch{1}.dtijobs.tracking.GTtrack.newfile.out.dir = {[options.root,options.patientname,filesep,'vat_resources',filesep]};
    matlabbatch{1}.dtijobs.tracking.GTtrack.newfile.out.fname = 'leadfield_FTR.mat';
    matlabbatch{1}.dtijobs.tracking.GTtrack.trackingarea.maskstruct.filenameMASK = {[options.root,options.patientname,filesep,'vat_resources',filesep,'trackmask_b0.mat']};
    matlabbatch{1}.dtijobs.tracking.GTtrack.trackingarea.maskstruct.roiid = 'trackmask_b0';
    matlabbatch{1}.dtijobs.tracking.GTtrack.parameters = 1;
    matlabbatch{1}.dtijobs.tracking.GTtrack.para_weight.custom_para_weight = 0.2;
    matlabbatch{1}.dtijobs.tracking.GTtrack.para_other.custom_para_other = [0.1 0.001 50 500000000 1 2 0.2 1];
    matlabbatch{1}.dtijobs.tracking.GTtrack.minlen = 1;
    matlabbatch{1}.dtijobs.tracking.GTtrack.maxlen = Inf;
    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear matlabbatch jobs
    disp('Done fibertracking.');
end

keyboard
%% generate headmodel from DTI fibers

%% generate leadfield:

%% old maedler code

[xx,yy,zz]=psphere(100);

try
    if isnan(stimparams)
        clear stimparams
    end
end
    radius=repmat(1.5,options.elspec.numel,1); % some default setting.
    %try % if stimparams are set.
    for con=1:length(stimparams(1,side).U)


        radius(con)=maedler12_eq3(stimparams(1,side).U(con),stimparams(1,side).Im(con));
        volume(con)=(4/3)*pi*radius(con)^3;

        VAT{con}=[xx*radius(con)+coords{side}(con,1);...
            yy*radius(con)+coords{side}(con,2);...
            zz*radius(con)+coords{side}(con,3)]';
    end
varargout{1}=VAT;
varargout{2}=radius;
varargout{3}=volume;


function r=maedler12_eq3(U,Im)
% This function radius of Volume of Activated Tissue for stimulation settings U and Ohm. See Maedler 2012 for details.
% Clinical measurements of DBS electrode impedance typically range from
% 500?1500 Ohm (Butson 2006).
r=0; %
if U %(U>0)

k1=-1.0473;
k3=0.2786;
k4=0.0009856;

r=-(k4*Im-sqrt(k4^2*Im^2  +   2*k1*k4*Im    +   k1^2 +   4*k3*U)   +   k1)...
    /...
    (2*k3);
end


function [x,y,z,avgr] = psphere(n,rad)
% [x,y,z,r] = psphere(N)
%
% Distributes N points "equally" about a unit sphere.
%
% N The number of points to distribute
% x,y,z Each is 1 x N vector
% r The smallest linear distance between two neighboring
% points. If the function is run several times for the
% same N, r should not change by more than the convergence
% criteria, which is +-0.01 on a unit sphere.
%
%
% Distributes N points about a unit sphere so that the straight line
% distance between neighboring points is roughly the same. The
% actual criteria for stopping is slightly different. The difference
% between a given point and every other point is stored in a matrix.
% The iterations stop once the maximum difference between any element
% in successive distance matrices is less than 0.01. An absolute
% criteria was chosen due to self-distances being 0, and programming
% around this for relative convergence seemed too much work for too
% little reward.
%
% The algorithm first generates N random points. Then a repulsive
% force vector, based on 1/r^2, is calculated for each point due to
% all of the other points. The resultant force vector for each point
% is normalized, and then each point is displaced a distance S = 1
% in the direction of the force. Any value of S between 0.0 and 1.
% 0 seems to work with most values between 0.2 and 1 taking an average
% of 20 to 30 iterations to converge. If S is too high, too much
% "energy" is being added to the system and it won't converge. If S is
% too low, the points won't be evenly distributed even though the
% convergence criteria is met. Generally speaking, the larger N is
% the larger S can be without causing instabilities. After
% displacement, the point is projected back down onto the unit sphere.
% When the system nears convergence, the displacement vector for a
% given point is nearly in the same direction as the radius vector for
% that point due to the points being equally distributed. A good check
% to make sure the code is working is to specify N = 4 and then check
% that the resulting points form a regular tetrahedron (or whatever
% it's called). How you would do this I don't know (check angles
% maybe). That's why I supplied the demo option so you could look at
% it in progress.
%
% Jason Bowman
% jbowman90@hotmail.com
% Last Revised June 2000
%
% You may freely distribute this code. I only ask that you give credit
% where credit is due.
%

%Since rand produces number from 0 to 1, subtract off -0.5 so that
%the points are centered about (0,0,0).
x = rand(1,n) - 0.5;
y = rand(1,n) - 0.5;
z = rand(1,n) - 0.5;

%Make the matrix R matrices for comparison.
rm_new = ones(n);
rm_old = zeros(n);

%Scale the coordinates so that their distance from the origin is 1.
r = sqrt(x.^2 + y.^2 + z.^2);

x = x./r;
y = y./r;
z = z./r;

not_done = 1;

s = 1;

%Turns off the divide by 0 warning
warning off

while not_done

    for i = 1:n

        %Calculate the i,j,k vectors for the direction of the repulsive forces.
        ii = x(i) - x;
        jj = y(i) - y;
        kk = z(i) - z;

        rm_new(i,:) = sqrt(ii.^2 + jj.^2 + kk.^2);

        ii = ii./rm_new(i,:);
        jj = jj./rm_new(i,:);
        kk = kk./rm_new(i,:);

        %Take care of the self terms.
        ii(i) = 0;
        jj(i) = 0;
        kk(i) = 0;

        %Use a 1/r^2 repulsive force, but add 0.01 to the denominator to
        %avoid a 0 * Inf below. The self term automatically disappears since
        %the ii,jj,kk vectors were set to zero for self terms.
        f = 1./(0.01 + rm_new(i,:).^2);

        %Sum the forces.
        fi = sum(f.*ii);
        fj = sum(f.*jj);
        fk = sum(f.*kk);

        %Find magnitude
        fn = sqrt(fi.^2 + fj.^2 + fk.^2);

        %Find the unit direction of repulsion.
        fi = fi/fn;
        fj = fj/fn;
        fk = fk/fn;

        %Step a distance s in the direciton of repulsion
        x(i) = x(i) + s.*fi;
        y(i) = y(i) + s.*fj;
        z(i) = z(i) + s.*fk;

        %Scale the coordinates back down to the unit sphere.
        r = sqrt(x(i).^2 + y(i).^2 + z(i).^2);

        x(i) = x(i)/r;
        y(i) = y(i)/r;
        z(i) = z(i)/r;

    end

    %Check convergence
    diff = abs(rm_new - rm_old);

    not_done = any(diff(:) > 0.01);

    rm_old = rm_new;
end %while

%Find the smallest distance between neighboring points. To do this
%exclude the self terms which are 0.
tmp = rm_new(:);

indices = find(tmp~=0);

avgr = min(tmp(indices));

%Turn back on the default warning state.
warning backtrace
