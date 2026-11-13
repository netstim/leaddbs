function tests = test_ea_getrawct2preniimat
% Raw-CT world coordinates must agree with actual ANTs image resampling.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
addpath(repo, fullfile(repo,'helpers'), fullfile(repo,'dependency'), ...
    fullfile(repo,'ext_libs'), fullfile(repo,'ext_libs','GetFullPath'), ...
    fullfile(repo,'ext_libs','iso2mesh'));
testCase.addTeardown(@() path(oldPath));
assumeTrue(testCase, exist('spm_vol','file') ~= 0, 'SPM must be on the MATLAB path.');
executable = ea_getExec(fullfile(repo,'ext_libs','ANTs','antsApplyTransforms'), ...
    escapePath = false);
assumeTrue(testCase,isfile(executable),'The bundled ANTs executable is required.');
testCase.TestData.ants = executable;
end

function setup(testCase)
folder = tempname;
mkdir(folder);
testCase.addTeardown(@() rmdir(folder,'s'));
testCase.TestData.folder = folder;
end

function testMatchingFormsAgreeWithAntsReslice(testCase)
f = fixture(testCase,'matching');
[tmat,ct] = ea_getrawct2preniimat(f.options);
verifyEqual(testCase,ct,f.raw);
verifyEqual(testCase,tmat,f.transform,'AbsTol',1e-10);
verifyAntsReslice(testCase,f,tmat);
end

function testGantryTiltSformDisagreementAgreesWithAntsReslice(testCase)
f = fixture(testCase,'sheared');
[tmat,ct] = ea_getrawct2preniimat(f.options);
verifyEqual(testCase,ct,f.raw);
verifyAntsReslice(testCase,f,tmat);
% These landmarks have known voxel positions. The expected physical points
% use the orthogonal scanner geometry selected by the actual ANTs reader.
landmarks = [2 7 15; 3 11 18; 2 8 13; 1 1 1];
actual = tmat * f.sform * landmarks;
expected = f.transform * f.qform * landmarks;
verifyEqual(testCase,actual,expected,'AbsTol',2e-4);
end

function testOrthogonalScannerSformIsNotReplacedByQform(testCase)
f = fixture(testCase,'scannerSform');
[tmat,~] = ea_getrawct2preniimat(f.options);
verifyAntsReslice(testCase,f,tmat);
landmarks = [2 7 15; 3 11 18; 2 8 13; 1 1 1];
verifyEqual(testCase,tmat * f.sform * landmarks, ...
    f.transform * f.sform * landmarks,'AbsTol',2e-4);
end

function testRoundedHalfTurnQuaternionUsesNiftiReferenceConvention(testCase)
f = fixture(testCase,'roundedQuaternion');
source = spm_vol(f.raw);
% Float32 rounding puts 1-sum(bcd.^2) below 1e-7. Reference NIfTI/ITK
% normalizes bcd and sets a=0; SPM instead leaves a small extra rotation.
verifyGreaterThan(testCase,max(abs(source.private.mat0-f.qform),[],'all'),1e-4);
[itkMat,spmMat] = ea_ants_image_affine(f.raw);
verifyEqual(testCase,itkMat,f.qform,'AbsTol',2e-6);
verifyEqual(testCase,spmMat,source.mat);
verifyEqual(testCase,itkMat([2 3],1),zeros(2,1),'AbsTol',1e-12);
verifyEqual(testCase,itkMat(1,[2 3]),zeros(1,2),'AbsTol',1e-12);
tmat = ea_getrawct2preniimat(f.options);
verifyAntsReslice(testCase,f,tmat);
end

function testInverseUndoesCorrectedForwardMapping(testCase)
f = fixture(testCase,'sheared');
[forward,ctForward] = ea_getrawct2preniimat(f.options,0);
[inverse,ctInverse] = ea_getrawct2preniimat(f.options,1);
verifyEqual(testCase,ctForward,f.raw);
verifyEqual(testCase,ctInverse,f.raw);
verifyEqual(testCase,inverse * forward,eye(4),'AbsTol',2e-10);
landmarks = [2 7 15; 3 11 18; 2 8 13; 1 1 1];
anchorPoints = f.transform * f.qform * landmarks;
verifyEqual(testCase,inverse * anchorPoints,f.sform * landmarks,'AbsTol',2e-4);
end

function testReslicedCtUsesIdentityInBothDirections(testCase)
f = fixture(testCase,'sheared');
f.options.prefs.reco.mancoruse = 'rpostop';
for inverse = [0 1]
    [tmat,ct] = ea_getrawct2preniimat(f.options,inverse);
    verifyEqual(testCase,tmat,eye(4));
    verifyEqual(testCase,ct,f.registered);
end
end

function testMissingCtMethodStillCorrectsExistingAntsTransform(testCase)
f = fixture(testCase,'sheared');
writeMethod(f.options.subj.coreg.log.method,false);
[tmat,ct] = ea_getrawct2preniimat(f.options);
verifyEqual(testCase,ct,f.raw);
verifyAntsReslice(testCase,f,tmat);
inverse = ea_getrawct2preniimat(f.options,1);
verifyEqual(testCase,inverse * tmat,eye(4),'AbsTol',2e-10);
end

function testMissingCtMethodAndTransformUsesReslicedCt(testCase)
f = fixture(testCase,'sheared');
writeMethod(f.options.subj.coreg.log.method,false);
delete([f.options.subj.coreg.transform.CT.forwardBaseName 'ants.mat']);
delete([f.options.subj.coreg.transform.CT.inverseBaseName 'ants.mat']);
for inverse = [0 1]
    [tmat,ct] = ea_getrawct2preniimat(f.options,inverse);
    verifyEqual(testCase,tmat,eye(4));
    verifyEqual(testCase,ct,f.registered);
end
end

function f = fixture(testCase,geometry)
folder = testCase.TestData.folder;
f.raw = fullfile(folder,'raw CT.nii');
f.reference = fullfile(folder,'anchor.nii');
f.registered = fullfile(folder,'registered CT.nii');
angle = 0.23;
rotation = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];
qzero = [rotation * diag([-0.7 1.1 2.3]), [13;-17;23]; 0 0 0 1];
if strcmp(geometry,'roundedQuaternion')
    coefficients = double(single([0 0.9930684566497803 0.11753722280263901]));
    normalized = coefficients / norm(coefficients);
    c = normalized(2);
    d = normalized(3);
    % This is the known a=0 reference quaternion, with negative qfac.
    rotation = [-1 0 0;0 c*c-d*d 2*c*d;0 2*c*d d*d-c*c];
    qzero(1:3,1:3) = rotation * diag([0.7 1.1 -2.3]);
end
szero = qzero;
switch geometry
    case {'sheared','roundedQuaternion'}
        % A CT gantry tilt can be represented in sform but not in qform.
        szero(1:3,3) = [0;0;2.3];
    case 'scannerSform'
        % An orthogonal scanner sform is valid and must take precedence.
        szero(1:3,4) = szero(1:3,4) + [3;-5;2];
end
oneToZero = [eye(3),-ones(3,1);0 0 0 1];
f.qform = qzero * oneToZero;
f.sform = szero * oneToZero;
[i,j,k] = ndgrid(1:18,1:21,1:15);
data = single(100 + i + 10*j + 30*k + 7*sin(i/3).*cos(j/4));
writeNifti(f.raw,data,f.sform,f.qform);
if strcmp(geometry,'roundedQuaternion')
    % Preserve the float32 quaternion coefficients exactly: a round trip
    % through SPM's mat0 encoder can change which threshold is crossed.
    fid = fopen(f.raw,'r+','ieee-le');
    assert(fid ~= -1,'Cannot patch fixture quaternion.');
    fseek(fid,256,'bof'); % NIfTI-1 quatern_b/c/d byte offsets 256/260/264.
    fwrite(fid,single(coefficients),'single');
    fclose(fid);
end
referenceZero = [diag([0.8 0.8 0.8]),[-10;-30;10];0 0 0 1];
referenceMat = referenceZero * oneToZero;
writeNifti(f.reference,zeros(80,80,80,'single'),referenceMat,referenceMat);
% A known nontrivial moving-to-fixed transform in RAS physical units.
angle = 0.13;
rotation = [cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
f.transform = [rotation * diag([1.03 0.96 1.02]),[2;-1;3];0 0 0 1];
forward = fullfile(folder,'from-CT_to-anchor_');
inverse = fullfile(folder,'from-anchor_to-CT_');
writeAntsTransform([forward 'ants.mat'],f.transform);
writeAntsTransform([inverse 'ants.mat'],inv(f.transform));
logfile = fullfile(folder,'coregmethod.json');
writeMethod(logfile,true);
f.options.prefs.reco.mancoruse = 'postop';
f.options.subj.preproc.anat.postop.CT = f.raw;
f.options.subj.coreg.anat.postop.CT = f.registered;
f.options.subj.coreg.log.method = logfile;
f.options.subj.coreg.transform.CT.forwardBaseName = forward;
f.options.subj.coreg.transform.CT.inverseBaseName = inverse;
end

function verifyAntsReslice(testCase,f,tmat)
transformFile = [f.options.subj.coreg.transform.CT.forwardBaseName 'ants.mat'];
cmd = sprintf('%s -d 3 -i %s -r %s -o %s -t %s -n Linear', ...
    ea_path_helper(testCase.TestData.ants),ea_path_helper(f.raw), ...
    ea_path_helper(f.reference),ea_path_helper(f.registered), ...
    ea_path_helper(transformFile));
[status,output] = system(cmd);
assertEqual(testCase,status,0,output);
source = spm_vol(f.raw);
registered = spm_vol(f.registered);
data = spm_read_vols(source);
actual = spm_read_vols(registered);
% Sample an interior grid of ANTs output voxels using the transform returned
% by Lead-DBS and SPM's own raw-image affine. This detects wrong q/s choice,
% wrong composition order, transform direction, and voxel-origin offsets.
[i,j,k] = ndgrid(2:2:78,2:2:78,2:2:78);
voxels = [i(:)';j(:)';k(:)';ones(1,numel(i))];
rawVoxels = source.mat \ (tmat \ (registered.mat * voxels));
values = actual(sub2ind(size(actual),i(:),j(:),k(:)));
inside = all(rawVoxels(1:3,:) > 1.1,1) & ...
    all(rawVoxels(1:3,:) < source.dim(:)-0.1,1) & values' > 0;
assertGreaterThan(testCase,nnz(inside),500);
sampled = interpn(data,rawVoxels(1,inside),rawVoxels(2,inside), ...
    rawVoxels(3,inside),'linear');
residual = sampled(:) - values(inside);
verifyLessThan(testCase,sqrt(mean(residual.^2)),0.01, ...
    'Raw CT sampled through Lead-DBS must reproduce the actual ANTs reslice.');
end

function writeNifti(filename,data,sform,qform)
volume = nifti;
volume.dat = file_array(filename,size(data),'FLOAT32-LE',352,1,0);
volume.mat = sform;
volume.mat0 = qform;
volume.mat_intent = 'Scanner';
volume.mat0_intent = 'Scanner';
create(volume);
volume.dat(:,:,:) = data;
end

function writeAntsTransform(filename,forwardRas)
rasToLps = diag([-1 -1 1 1]);
pullLps = rasToLps / forwardRas * rasToLps;
fixed = [5;-3;2];
translation = pullLps(1:3,4) - fixed + pullLps(1:3,1:3)*fixed;
AffineTransform_double_3_3 = [reshape(pullLps(1:3,1:3)',9,1);translation];
% ITK uses MATLAB v4 with full transform-type names. MATLAB's save -v4
% refuses names longer than 19 characters, although its reader and ITK
% support them, so write the two plain real-valued v4 arrays directly.
fid = fopen(filename,'w','ieee-le');
assert(fid ~= -1,'Cannot write fixture ANTs transform.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
writeMat4Array(fid,'AffineTransform_double_3_3',AffineTransform_double_3_3);
writeMat4Array(fid,'fixed',fixed);
end

function writeMat4Array(fid,name,values)
fwrite(fid,int32([0 size(values,1) size(values,2) 0 numel(name)+1]),'int32');
fwrite(fid,[name char(0)],'char');
fwrite(fid,values,'double');
end

function writeMethod(filename,withCt)
fid = fopen(filename,'w');
assert(fid ~= -1,'Cannot write fixture method log.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
if withCt
    fprintf(fid,'{"method":{"CT":"ANTs (Avants 2008)"}}');
else
    fprintf(fid,'{"method":{"T2w":"SPM (Friston 2007)"}}');
end
end
