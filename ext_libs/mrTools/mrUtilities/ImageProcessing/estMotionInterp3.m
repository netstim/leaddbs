function M = estMotionInterp3(tseries1,tseries2,frame1,frame2,numIters,Minitial,sliceTimes,rotFlag,robustFlag,phaseFlag,crop,CB,SC)
%
% function M = estMotionInterp3(tseries1,tseries2,frame1,frame2,...
%              [numIters],[Minitial],[sliceTimes],[rotFlag],[robustFlag],[phaseFlag],[crop])
%
% tseries1 and tseries2 are time series of volumes (4d arrays)
%     If either tseries has only 1 frame then that frame is used with no
%     slice time correction.
% frame1 and frame2 are corresponding frames from each of the tseries that
%     are being compared.
% numIters is number of iterations to run
% Minitial is initial guess for M. Default is 3x3 identity matrix.
% sliceTimes: vector of length equal to the number of slices. Default is a
%     vector of zeros.
% crop specifies border size to crop/ignore  around all sides of the volume. 
%     Should be of the form [ymin xmin zmin; ymax xmax zmax]
%     Default crops 2 pixel border: [2 2 2; (size(vol1) - [2 2 2])].
% rotFlag: passed to estMotion3 (if True, contrains rigid body motion).
%     Default: 1
% robustFlag is passed to estMotion3 (if True, uses robust
%     M-estimator). Default: 0.
% phaseFlag: If True, treats the input volumes as containing the phase of
%     complex-values (on the unit circle in the complex plane) and computes
%     the motion estimates on the complex values. This is useful if the
%     image intensities in the volumes are from a periodic domain such that
%     the "brightest" and "darkest" intensities are actually very similar
%     to each other. Default: 0.
%
% M is 4x4 translation+rotation or affine transform matrix: X' = M X
% where X=(x,y,1) is starting position in homogeneous coords
% and X'=(x',y'',1) is ending position
%
% where X=(x,y,z,1) is starting position in homogeneous coords
% and X'=(x',y',z',1) is ending position
%
% Each iteration performs slice time correction while warping the volumes
% according to the previous estimate, and then estimates the residual
% motion.
%
% DJH 7/2007, modified from estMotionIter3

% default values
tseries1Size = size(tseries1);
if ieNotDefined('numIters')
    numIters=3;
end
if ieNotDefined('sliceTimes')
  sliceTimes = zeros(1,tseries1Size(3));
end
if ieNotDefined('Minitial')
    M=eye(4);
else
    M=Minitial;
end
if ieNotDefined('rotFlag')
    rotFlag = 1;
end
if ieNotDefined('robustFlag')
    robustFlag = 0;
end
if ~exist('phaseFlag','var')
    phaseFlag = 0;
end
if ieNotDefined('crop')
    crop = [2 2 2; (tseries1Size(1:3) - [2 2 2])];
end
if ieNotDefined('CB')
    CB = [];
end
if ieNotDefined('SC')
    SC = [];
end

nframes1 = size(tseries1,4);
nframes2 = size(tseries2,4);

for iter=1:numIters
  Mhalf2=real(sqrtm(M));
  Mhalf1=real(sqrtm(inv(M)));
  if (nframes1 == 1)
    volWarp1=warpAffine3(tseries1(:,:,:,1),Mhalf1);
  else
    volWarp1=warpAffineInterp3(tseries1,frame1,Mhalf1,sliceTimes);
  end
  if (nframes2 == 1)
    volWarp2=warpAffine3(tseries2(:,:,:,1),Mhalf2); %JB: I think this is supposed to be tseries2 and not tseries1 (although this case probably never happens)
  else
    volWarp2=warpAffineInterp3(tseries2,frame2,Mhalf2,sliceTimes);
  end
  deltaM=estMotion3(volWarp1,volWarp2,rotFlag,robustFlag,phaseFlag,crop,CB,SC);
  M=deltaM*M;
end

return;


%%%%%%%%%
% Debug %
%%%%%%%%%

% Simple translation with no slice time correction, random image
% intensities
filter = [0.03504 0.24878 0.43234 0.24878 0.03504];
in = convXYsep(convZ(rand(30,40,20),filter),filter,filter);
A= [1 0 0 .2;
    0 1 0 .3;
    0 0 1 .4;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
tseries=zeros([size(in),2]);
tseries=zeros([size(in),2]);
tseries(:,:,:,1) = vol1;
tseries(:,:,:,2) = vol2;
estMotionInterp3(tseries,tseries,1,2,3)
estMotionInterp3(vol1,tseries,1,2,3)
estMotionInterp3(vol1,vol2,1,1,3)


% zone plate volumes
% tseries1 has slice time offsets
% tseries2 has no slice time correction
Nx=32;
Ny=34;
Nz=20;
Nt=60;
tseries1 = zeros(Ny,Nx,Nz,Nt);
tseries2 = zeros(Ny,Nx,Nz,Nt);
for xs = 1:Nx
  for ys = 1:Ny
    for zs = 1:Nz
      for ts = 1:Nt
        shift1 = ((ts-1)*Nz + zs)/(Nz*Nt) * 4* Nx;
        shift2 = ((ts-1)*Nz + 1)/(Nz*Nt) * 4* Nx;
        x1 = 2*pi*(4/Nx)*(xs-shift1);
        x2 = 2*pi*(4/Nx)*(xs-shift2);
        y = 2*pi*(4/Ny)*ys;
        z = 2*pi*(4/Ny)*zs;
        tseries1(ys,xs,zs,ts) = sin(x1)+sin(y)+sin(z);
        tseries2(ys,xs,zs,ts) = sin(x2)+sin(y)+sin(z);
      end
    end
  end
end

% views slices through the volumes
zs = 1;
for ts = 1:Nt
  figure(1);
  imshow(tseries1(:,:,10,ts),[-3 3],'InitialMagnification',200);
  figure(2);
  plot([squeeze(tseries1(16,:,zs,ts)),squeeze(tseries2(16,:,zs,ts))]);
  ylim([-3 3]);
  pause(0.1)
end

% Estimate motion
Minitial = eye(4);
Niters = 5;
sliceTimes1 = [0:Nz-1]/Nz;
sliceTimes2 = zeros(1,Nz);
frame1 = 2;
frame2 = 3;
expectedShift = (frame2-frame1) * (Nz)/(Nz*Nt) * 4 * Nx

% Estimate with correct slice times
estMotionInterp3(tseries1,tseries1,frame1,frame2,Niters,Minitial,sliceTimes1,1)
estMotionInterp3(tseries2,tseries2,frame1,frame2,Niters,Minitial,sliceTimes2,1)

% Estimate with incorrect slice times
estMotionInterp3(tseries1,tseries1,frame1,frame2,Niters,Minitial,sliceTimes2,1)
estMotionInterp3(tseries2,tseries2,frame1,frame2,Niters,Minitial,sliceTimes1,1)






