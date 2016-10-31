% fid2xform.m
%
%        $Id: fid2xform.m 2855 2013-09-05 06:51:33Z justin $ 
%      usage: [xform info] = fid2xform(fidname,<verbose>)
%         by: justin gardner
%       date: 05/05/09
%    purpose: Convert fields out of a named procpar into a rotation matrix
%             suitable for the nifti header. (fidname is either the name of 
%             a fid, or the structure returned by readprocpar). info is 
%             a sturcture that contains info about the scan like pixel dims
%             sense factors etc.
%
function [xform info] = fid2xform(procpar,verbose,varargin)

xform = [];info = [];

% check arguments
if nargin < 1
  help fid2xform
  return
end

% get extra arguments
movepro=[];movepss=[];
getArgs(varargin,{'movepro=0','movepss=0','fixConsoleUpgradeBug=1'});

if ieNotDefined('verbose'),verbose = 0;end

% get the procpar
if isstr(procpar)
  info.fidname = procpar;
  procpar = readprocpar(procpar);
  if isempty(procpar),return,end;
elseif ~isstruct(procpar)
  help fid2xform;
  return
end

% get the slice order
if isfield(procpar,'pss')
  [sortedpss sliceOrder] = sort(procpar.pss);
else
  sliceOrder = [];
end

% move pro if called for
if movepro ~= 0
  procpar.pro = procpar.pro + movepro;
end

% get the number of navechoes
if isfield(procpar,'navechoes')
  navechoes = procpar.navechoes;
else
  navechoes = 0;
end

% check to see if this is an epi
info.isepi = 0;
if (isfield(procpar,'petable'))
  token = procpar.petable{1};
  % the petable name should be something like
  % "epi132alt8k". If so, it is an epi.
  if (strncmp(token,'epi',3)) 
    info.isepi = 1;
  % some files may not be named according to convention
  % in that case we check for the cntr field. If it is there
  % then this is an epi
  elseif isfield(procpar,'cntr') && ~isequal(procpar.cntr,0)
    info.isepi = 1;
  end
end


% get date and time of start and finish
[info.startDatevec info.startDatestr] = getDateFromVarianField(procpar.time_run);
[info.endDatevec info.endDatestr] = getDateFromVarianField(procpar.time_complete);
% get elapsed time
info.elapsedSecs = etime(info.endDatevec,info.startDatevec);
elapsedMin = floor(info.elapsedSecs/60);
info.elapsedTimeStr = sprintf('%i min %i sec',elapsedMin,info.elapsedSecs-elapsedMin*60);

% set the console in the info for easy reference
info.console = procpar.console{1};

% get the intlv setting
if isfield(procpar,'intlv')
  if iscell(procpar.intlv)
    info.intlv = procpar.intlv{1};
  else
    info.intlv = procpar.intlv;
  end
else
  if isfield(procpar,'transhim')
    info.intlv = 'y';
    procpar.numshots = procpar.nv;
  else
    info.intlv = 'n';
  end
end

% check for an error in intlv setting, which happens when intlv is set to 'y' and
% this is a 3D mprage (which shouldn't happen)
if strcmp(info.intlv,'y') && isequal(procpar.nD,3)
  disp(sprintf('(fid2xform) !!! Intlv is set to y even though this is a 3D mprage. This is likely a mistake. Setting intlv to n. !!!'));
  info.intlv = 'n';
end

% get the ilts setting
if isfield(procpar,'ilts')
  if iscell(procpar.ilts)
    info.ilts = procpar.ilts{1};
  else
    info.ilts = procpar.ilts;
  end
else
  info.ilts = [];
end

% if there is no ni field then just return
if ~isfield(procpar,'ni')
  disp(sprintf('(fid2xform) ni is not set. This may be a non-image file'))
  info.dim = nan(1,4);
  return
end

% get the dimensions of the scan
% (procpar.ni is lines of k-space, procpar.nv is number of lines collected including navigator echoes)
% used to use procpar.nv and correct for navechoes, but seems more sensible to just use procpar.ni)
dim = [procpar.np/2 procpar.ni length(procpar.pss)];

if (procpar.ni == 1) && (procpar.nf > 1)
  if verbose>0
    disp(sprintf('(fid2xform) Fid file looks like it is compressed. Using nf for the 2nd dim'));
  end
  info.compressedFid = true;
  dim(2) = procpar.nf;
  % if interleaving then the second dimension has been multiplied by the
  % number of slices, so undo that
  if isequal(info.intlv,'y')
    dim(2) = dim(2)/dim(3);
  end
else
  info.compressedFid = false;
end

% remove navigator echoes from k-space 
%dim(2) = dim(2) - navechoes;

% check to see if this is a sense reconstruction and what the sense acceleration factor is
if procpar.accfactor > 1
  % fix dimensions
  dim(2) = dim(2)*procpar.accfactor;
  % print message
  if verbose>0
    disp(sprintf('(fid2xform) Found a sense factor of %i. Dims are now [%i %i]',procpar.accfactor,dim(1),dim(2)));
  end
end

% check for weird psi theta or phi
if length(procpar.psi) > 1
  disp(sprintf('(fid2xform) psi is set to array %s. resetting to %f',mlrnum2str(procpar.psi),procpar.psi(1)));
  procpar.psi = procpar.psi(1);
end
if length(procpar.theta) > 1
  disp(sprintf('(fid2xform) theta is set to array %s. resetting to %f',mlrnum2str(procpar.theta),procpar.theta(end)));
  procpar.theta = procpar.theta(end);
end
if length(procpar.phi) > 1
  disp(sprintf('(fid2xform) phi is set to array %s. resetting to %f',mlrnum2str(procpar.phi),procpar.phi(2)));
  procpar.phi = procpar.phi(2);
end

% make the rotation matrix from the procpar angles
rotmat = euler2rotmatrix(procpar.psi,-procpar.theta,-procpar.phi);

info.processed = 1;
if dim(2) == 0
 info.processed = 0;
 if verbose > 0
   disp(sprintf('(fid2xform) Unprocessed fid directory. Has dim(2)=0'));
 end
 voxsize = [nan nan nan 1];
 voxspacing = [nan nan nan 1];
else
  % get voxel sizes and put in diagonal to multiply against rotmat
  % so that we get the proper spacing
  voxsize = [10*procpar.lro/dim(1) 10*procpar.lpe/dim(2) procpar.thk 1];

  % vox spacing can be *different* from voxsize, if you skip in your pss
  voxspacing = [10*procpar.lro/dim(1) 10*procpar.lpe/dim(2) 10*median(diff(sort(procpar.pss))) 1];
  
  % if we have only one slice, then get vox spacing from thk (since there isn't really any vox spacing
  if length(procpar.pss) == 1
    voxspacing(3) = procpar.thk;
  end
end

% check for 3d acquisition
info.acq3d = 0;
if procpar.nv2 > 1
  if verbose>0,disp(sprintf('(fid2xform) 3D acquisition'));end
  % since the 3rd dimension is taken as a single slice with multiple
  % phase encodes, we have to get the voxel size and dimensions differently
  voxsize(3) = 10*procpar.lpe2/procpar.nv2;
  voxspacing(3) = 10*procpar.lpe2/procpar.nv2;
  dim(3) = procpar.nv2;
  % keep in structure that this is a 3d acquisition
  info.acq3d = 1;
  % check to see if it has been processed or not (i.e. pss should be of
  % correct length
  if length(procpar.pss) ~= procpar.nv2
    info.processed = 0;
  end
end


% checking for a small bug in the 3d images recosntruction introduced
% by fft3rd which shifted the data by interger amounts in the 
% "slice" dimension, rather than the exact amount specified by
% the pss. Correcting for that here if the data were processed
% with the old version of the program
if info.acq3d & info.processed
  if ~isfield(procpar,'fftw3dexe_processed') 
    info.processed = 0;
  else
    if procpar.fftw3dexe_processed == 0
      % get number of slices
      nSlices = length(procpar.pss);
      if iseven(nSlices)
	% find the midpoint of the pss array, this is where the center of the 3d
	% slab was set to.
	mid_pss = mean(procpar.pss);
	% get the distance in between each slice
	slicediff = voxspacing(3)/10;
	% now reset the pss, to what it actually should have been
%	procpar.pss = sort(mid_pss-slicediff*nSlices/2+slicediff:slicediff:mid_pss+slicediff*nSlices/2,2,'descend');
        procpar.pss = sort(mid_pss-slicediff*nSlices/2:slicediff:mid_pss+slicediff*nSlices/2-slicediff/2,2,'descend');
	% The following works as well, but the original pss array is rounded to significant digits, so 
	% might as well recompute like above and get all significant digits
%	procpar.pss = procpar.pss-slicediff/2; 
	if verbose >= 0
	  disp(sprintf('(fid2xform) Adjusting pss array to account for non-integer shift',length(procpar.pss)));
	end
      end
    end
  end
end

% set field that says whether fft2rdexe_processing happend or not
if isfield(procpar,'fftw3dexe_processed')
  info.fftw3dexe_processed = 1;
else
  info.fftw3dexe_processed = 0;
end

% move pss if called for
if movepss ~= 0
  % apply the compensation only if this is a 3D file
  if ~info.acq3d || info.fftw3dexe_processed
    disp(sprintf('(fid2xform) !!! Unable to move pss for 2D data. Ignoring desired pss shift of: %f !!!',movepss));
  else
    procpar.pss = procpar.pss - movepss*10;
  end
end

% check to see if this is an uncompressedFid and 3D in which
% case the pss just contains the slice center, so we need to 
% adjust it here
if info.compressedFid && info.acq3d && (length(procpar.pss) == 1)
  % set to automatically movepss since origin is center of magnet not center of slices
  info.movepss = -procpar.pss;
  info.movepss = 0;
  if verbose > 0,disp(sprintf('(fid2xform) Compressed 3D fid, pss of center of slab: %f',procpar.pss));end
  % compute location of first and last slice
  firstSlice = procpar.pss - voxspacing(3)*(procpar.nv2-1)/2;
  lastSlice = procpar.pss + voxspacing(3)*(procpar.nv2-1)/2;
  % now make array
  procpar.pss = (firstSlice:voxspacing(3):lastSlice)/10;
else
  info.movepss = 0;
end

% count number of receivers
if isfield(procpar,'rcvrs')
  % count the number of receivers that have been turned on
  info.numReceivers = length(strfind(procpar.rcvrs{1},'y'));
else
  info.numReceivers = 1;
end


% Now get the offset in mm from the center of the bore that the center of the
% volume is. We can not change the phase encode center. Note that dimensions
% are given in cm, so we must convert to mm.
% 
% Note here about 3D images. The pss is in "reverse" order in the sense that it
% goes from positive numbers to negative numbers. This is fixed in fid2nifti
% since fid2nifti sorts the pss and reorders the slices accordingly. This reordering fixes
% interleaved acquisition as well. Knowing that this will be the case,
% means that the offset to the first slice should always be the min(procpar.pss)
offset = 10*[-procpar.pro 0 min(procpar.pss)];

% make into a translation matrix
offset = [eye(3) offset';0 0 0 1];

% get the distance to the image origin (ie voxel 0,0,0) in number of voxels
% (i.e. the image dimensions divided by 2 as we assume that the offset specifies
% where the center of the volume is - (the third dimension is set by the
% pss, so we ignore it)
% note that I would have thought that we need to subtract 1 so that
% we get to the center of the voxel, but not subtracting 1 empirical
% seems to be correct - actually this may be because of the way the fft works though
%originOffset = -(dim-1)/2;
originOffset = -dim/2;
originOffset(3) = 0;
originOffset = [eye(3) originOffset';0 0 0 1];

% this swaps the dimensions to the coordinate frame that Nifti is expecting.
swapDim =[0 0 1 0;1 0 0 0;0 1 0 0;0 0 0 1];

% Another final fix. This gets left/right correct
if strcmp(lower(info.console),'inova')  
  swapDim2 =[0 0 -1 0;0 1 0 0;1 0 0 0;0 0 0 1];
else
  swapDim2 =[0 0 -1 0;0 1 0 0;-1 0 0 0;0 0 0 1];
end

% epi images appear to nead a flip in X and Y
if info.isepi
  % strange 1 voxel shifts needed after console update
  shiftBlech = eye(4);
  % before about 2011/09/01 the readout was flipped too
  %epiFlip = [-1 0 0 1;0 -1 0 1;0 0 1 0;0 0 0 1];
  if strcmp(lower(info.console),'inova')  
    epiFlip = eye(4);
  else
    % new console needs a phase encode flip
    if isfield(procpar,'pe_dir')
      pe_dir=procpar.pe_dir{1};
      if (pe_dir(1) == 'l')
        epiFlip = eye(4);
        % we also seem need these strange 1 voxel shifts to align with mprage data, blech, blech, blech.
        if fixConsoleUpgradeBug
          shiftBlech(1,4) = -1;
          shiftBlech(2,4) = -1;
        end
      else
       epiFlip = [1 0 0 1;0 -1 0 1;0 0 1 0;0 0 0 1];
        % we also seem need these strange 1 voxel shifts to align with mprage data, blech, blech, blech.
        if fixConsoleUpgradeBug
          shiftBlech(1,4) = -1;
          shiftBlech(2,4) = 1;
        end
      end
    else
      epiFlip = [1 0 0 1;0 -1 0 1;0 0 1 0;0 0 0 1];
        % we also seem need these strange 1 voxel shifts to align with mprage data, blech, blech, blech.
        if fixConsoleUpgradeBug
          shiftBlech(1,4) = -1;
          shiftBlech(2,4) = 1;
        end
    end
  end
  % if ppe is set then that is the shift in the phase encode direction (this is ignored by mprage)
  offset(2,4) = procpar.ppe*10;
else
  epiFlip = eye(4);
  shiftBlech = eye(4);
end

% now create the final shifted rotation matrix
xform = swapDim2*rotmat*swapDim*offset*diag(voxspacing)*epiFlip*shiftBlech*originOffset;

% testing rotmat
%rotmatpsi = euler2rotmatrix(procpar.psi,0,0);
%rotmattheta = euler2rotmatrix(0,procpar.theta,0);
%rotmatphi = euler2rotmatrix(0,0,-procpar.phi);
%xform = rotmatpsi*rotmattheta*rotmatphi*sliceOffset*swapDim*offset*diag(voxsize)*originOffset;

% round-off to zero
xform((xform < 1e-10) & (xform > - 1e-10)) = 0;

% verbose display only
if verbose > 0
  % display some info.
  disp(sprintf('(fid2xform) psi=%0.2f phi=%0.2f theta=%0.2f',procpar.psi,procpar.phi,procpar.theta));
  disp(sprintf('(fid2xform) Scan dims=[%i %i %i]',dim(1),dim(2),dim(3)));
  disp(sprintf('(fid2xform) Voxel size: [%0.2f %0.2f %0.2f]',voxsize(1),voxsize(2),voxsize(3)));
  disp(sprintf('(fid2xform) Voxel spacing: [%0.2f %0.2f %0.2f]',voxspacing(1),voxspacing(2),voxspacing(3)));
  disp(sprintf('(fid2xform) First slice offset: [%0.2f %0.2f %0.2f]',offset(1,4),offset(2,4),offset(3,4)));
  disp(sprintf('(fid2xform) pss = %s',mlrnum2str(procpar.pss)));
  disp(sprintf('(fid2xform) length of pss = %i',length(procpar.pss)));
  disp(sprintf('(fid2xform) offset to origin: [%0.2f %0.2f %0.2f]',originOffset(1),originOffset(2),originOffset(3)));
end

% work out what the tr is, this is just for setting in the info output field and
% is used by fid2niftihdr
% get the volume TR (framePeriod) for EPI 
tr = procpar.tr;
% if we run mutliple shots, volume TR = slice TR * shots 
if isfield(procpar,'navechoes')
  tr = tr * procpar.numshots/procpar.accfactor;
end
% if we ran without interleaving slices then
%  volume TR = slice TR * shots * slice number
if isfield(procpar,'intlv') && strcmp(procpar.intlv, 'n')
  tr = tr * length(procpar.pss);
end

% check for 2d anatomy, to tell getfid that this has the slices and receivers mixed up
if ~info.acq3d && ~info.isepi
  if strcmp(lower(info.console),'inova')  
    info.receiversAndSlicesSwapped = 1;
    if verbose>0,disp(sprintf('(fid2xform) Receviers and slices are swapped'));,end
  else
    % looks like on the new console we don't have to do this anymore?
    info.receiversAndSlicesSwapped = 0;
  end
else 
  info.receiversAndSlicesSwapped = 0;
end
  

% pack up some useful information that we have learned about the scan
% into an output structure that can be used by other programs.
info.dim = dim;
info.voxsize = voxsize;
info.voxspacing = voxspacing;
info.offset = offset;
info.originOffset = originOffset;
info.pss = procpar.pss;
info.psi = procpar.psi;
info.phi = procpar.phi;
info.theta = procpar.theta;
info.accFactor = procpar.accfactor;
info.nRefVolumes = 0;
if strcmp(lower(info.console),'inova') && isfield(procpar,'cntr')
  % compute the length of the scan minus the number of steady state reference volumes
  % are taken. This is computed by looking for the first volume of the scan that
  % has a trigger in it. 
  info.nRefVolumes = first(find(procpar.cntr))-1;
  info.dim(4) = length(procpar.cntr)-info.nRefVolumes;
elseif strcmp(lower(info.console),'vnmrs') && isfield(procpar,'image')
  % compute the length of the scan minus the number of steady state reference volumes
  % are taken. This is computed by looking for the first volume of the scan that
  % has a trigger in it. 
  info.nRefVolumes = first(find(procpar.image))-1;
  info.dim(4) = length(procpar.image)-info.nRefVolumes;
else
  info.dim(4) = 1;
end
info.tr = tr;
info.sliceOrder = sliceOrder;

% keep procpar
info.procpar = procpar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getDateFromVarianField    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outDateVec outDateStr] = getDateFromVarianField(f)

outDateVec = [];outDateStr = '';
if isempty(f{1}),return,end

year = str2num(f{1}(1:4));
month = str2num(f{1}(5:6));
day = str2num(f{1}(7:8));

hour = str2num(f{1}(10:11));
min = str2num(f{1}(12:13));
sec = str2num(f{1}(14:15));

outDateVec = datevec(sprintf('%04i/%02i/%02i %02i:%02i:%02i',year,month,day,hour,min,sec));
outDateStr = datestr(outDateVec);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   euler2rotmatrix   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function r = euler2rotmatrix(a,b,g)

% convert to radians
a = pi*a/180;b = pi*b/180;g = pi*g/180;

% get cos and sin of angles
ca = cos(a);sa = sin(a);
cb = cos(b);sb = sin(b);
cg = cos(g);sg = sin(g);

% convert each rotation into a rotation matrix
%arot1 = [ca  sa 0;-sa ca 0;0   0  1];
%arot2 = [ca  0 sa;0 1  0;-sa 0 ca];
arot3 = [1  0 0;0 ca  sa;0  -sa ca];
brot1 = [cb sb 0;-sb cb  0;0 0 1];
%brot2 = [cb 0 sb;0 1  0;-sb 0 cb];
%brot3 = [1 0   0; 0 cb  sb;0  -sb cb];
%grot1 = [cg  sg 0;-sg cg 0;0   0  1];
%grot2 = [cg  0  sg;0   1   0;-sg  0 cg];
grot3 = [1 0   0;0 cg  sg;0  -sg cg];

% composite the rotations together
r = arot3*brot1*grot3;

% make into a homogenized xform
r = [[r;0 0 0] [0 0 0 1]'];

