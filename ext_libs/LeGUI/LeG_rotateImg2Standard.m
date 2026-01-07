function [NiiFileOut,FlipB] = LeG_rotateImg2Standard(NiiFile,varargin)
%Performs 90 deg rotations to place image in standard space with dim1 (x ->
%saggital plane), dim2 (y -> coronal plane), dim3 (z -> axial plane). Finds
%matrix to transform nifti affine matrix to identity matrix with 1st dim
%mirrored ([-1 0 0; 0 1 0; 0 0 1]). Saves with "sd" prefix.
%
%Tyler Davis
%20220310

imgInfo = spm_vol(NiiFile);
img = spm_read_vols(imgInfo);
imgSize = size(img);
imgCenter = round(imgSize/2);
imgScale = sqrt(sum(imgInfo.mat(1:3,1:3).^2));
% imgCenterMM = imgCenter.*imgScale;

% imgInfoNew = imgInfo;

[path,file,ext] = fileparts(imgInfo.fname);
if isempty(regexp(file,'^sd','once'))
    imgInfo.fname = fullfile(path,['sd',file,ext]);
%     imgInfoNew.fname = fullfile(path,['nw',file,ext]);
end

%%%%%%%%%%%%%%%% Pre-rotation %%%%%%%%%%%%%%%%
% x = imgInfo.mat(1:3,1); [~,xi] = max(abs(x)); xs = sign(x(xi));
% y = imgInfo.mat(1:3,2); [~,yi] = max(abs(y)); ys = sign(y(yi));
% z = imgInfo.mat(1:3,3); [~,zi] = max(abs(z)); zs = sign(z(zi));
% 
% nmat = zeros(4); %new (desired) matrix
% nmat(xi,1) = imgScale(1)*xs;
% nmat(yi,2) = imgScale(2)*ys;
% nmat(zi,3) = imgScale(3)*zs;
% nmat(1,4) = imgCenterMM(1)*sign(imgInfo.mat(1,4));
% nmat(2,4) = imgCenterMM(2)*sign(imgInfo.mat(2,4));
% nmat(3,4) = imgCenterMM(3)*sign(imgInfo.mat(3,4));
% nmat(4,4) = 1;
% 
% % rmat = nmat\imgInfo.mat; %mapping from imgInfo.mat into nmat (vox2vox)
% 
% flags.mask = 1;
% flags.interp = 4; %4th order bspline
% flags.which = [1,0];
% flags.wrap = [0,0,0];
% flags.prefix = '';
% 
% imgInfo.mat(1:3,4) = nmat(1:3,4); %might want the translations to be the same and new matrix
% imgInfoNew.mat = nmat; %write desired matrix to header
% 
% disp([imgInfo.mat,nan(4,1),nmat]);
% mmvox = sqrt(sum((imgInfo.mat*ones(4,1)-(imgInfo.mat*(diag(ones(1,4))+1))).^2)); %calculating distance in mm between adjacent voxels along each dimension
% mmvoxnew = sqrt(sum((nmat*ones(4,1)-(nmat*(diag(ones(1,4))+1))).^2));
% disp([mmvox(1:3),nan,mmvoxnew(1:3)])
% 
% spm_write_vol(imgInfo,img); %create the nifti file (sdMR.nii) with the original matrix
% spm_write_vol(imgInfoNew,img); %create the nifti file (orMR.nii) with the new matrix
% 
% P(1) = {[imgInfoNew.fname,',1']}; %contains the desired orientation
% P(2) = {[imgInfo.fname,',1']}; %contains the sd prefix and the original orientation (this one is resliced)
% 
% fHProg = findobj('tag','Interactive');
% if isempty(fHProg)
%     fHProg = figure('tag','Interactive'); %progress figure
% end
% spm_reslice(P,flags);
% 
% imgInfo = spm_vol(imgInfo.fname);
% img = spm_read_vols(imgInfo);
% 
% nSize = size(img);
% nScale = sqrt(sum(imgInfo.mat(1:3,1:3).^2));
% if any(nSize~=imgSize)||any(round(nScale,4)~=round(imgScale,4))||any(any(round(imgInfo.mat,4)~=round(nmat,4)))
%     msgbox('Pre-rotation failed in LeG_rotateImg2Standard');
% end
% delete(imgInfoNew.fname);  
% 
% close(fHProg); pause(0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Performing 90 deg rotations to put in standard space
eyemat = eye(4); 
eyemat(1,1) = -1;
rmat = imgInfo.mat\eyemat; %imgInfo.mat*rmat = eyemat (this might be reversed - eyemat\imgInfo.mat is voxvox transform from imgInfo to eyemat)
b = spm_imatrix(rmat);
rvec = b(4:6);

imgR = img;
nrot = round(rvec/(pi/2)); %number of rotations
rvec = nrot*(pi/2); %rounded to nearest 90 deg rotation
for m=1:3 %spm order (shear, scale, rotation, translation)
    pmat = [setdiff(1:3,m),m];
    if nrot(m)~=0
        imgR = permute(imgR,pmat);
        imgR = rot90(imgR,nrot(m));
        imgR = ipermute(imgR,pmat);
    end
end

b0 = zeros(1,12);
b0(4:6) = rvec;
b0(7:9) = 1;
rmat0 = spm_matrix(b0); %this is now the rotation matrix rounded to nearest 90 deg

trans = abs((imgCenter.*imgScale)*rmat0(1:3,1:3)); trans(2:3) = -trans(2:3);
scl = abs(imgScale*rmat0(1:3,1:3)); scl(1) = -scl(1);

eyemat = eye(4);
eyemat(1:3,1:3) = diag(scl);
eyemat(1:3,4) = trans;

%%%%%%%%%%%%%%%%%%%%% check if a flip is needed %%%%%%%%%%%%%%%%%%%%%
imat = imgInfo.mat*rmat0;
x = imat(1:3,1); [~,xi] = max(abs(x)); xs = sign(x(xi));
FlipB = false;
if xs==1
    FlipB = true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgInfo.mat = eyemat;
imgInfo.dim = size(imgR);
spm_write_vol(imgInfo,imgR);

NiiFileOut = imgInfo.fname;