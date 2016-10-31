% getr2timecourse.m
%
%        $Id$
%      usage: d = getr2timecourse(timecourses,nhdr,hdrlen,scm,<tr>)
%         by: justin gardner
%       date: 05/25/07
%    purpose: function that will use getr2 to do the event related
%             analysis. TR argument is only necessary to return
%             also the time vector. The timecourses argument is
%             a kxn matrix with the length n timecourses for k
%             voxels. nhdr (number of hemodynamic responses) hdrlen
%             is the lenght of the hemodynamic response, and scm
%             is the stimulus convolution matrix. usually these
%             will be taken from the "d" variable returned from
%             the eventRelated code.
% 
%
function outd = getr2timecourse(timecourses,nhdr,hdrlen,scm,tr,verbose)

% check arguments
if ~any(nargin == [4 5 6])
  outd = [];
  help getr2timecourse
  return
end

if ieNotDefined('verbose'),verbose = 1;end

% make sure this has non-zero mean (otherwise event-related code
% which converts to percent signal change will divide by a value
% close to zero which causes the responses to to have strange magnitudes)
if max(abs(mean(timecourses,2))) < 1e-10
  if verbose
    disp(sprintf('(getr2timecourse) Adding 1 to timecourses because timecourses apparently have 0 mean'));
  end
  timecourses = timecourses+1;
end

% make sure every dimension has more than one timecourse
% even if we have to repeat.
for i = 1:2
  for j = 1:2
    if size(timecourses,1) == 1
      d.data(i,j,1,:) = timecourses;
      d.data(i,j,2,:) = timecourses;
    else
      d.data(i,j,:,:) = timecourses;
    end
  end
end

% set fields
d.scm = scm;
d.hdrlen = hdrlen;
d.nhdr = nhdr;
d.dim = size(d.data);
d.volumes = 1:d.dim(4);

% compute the event related analysis
outd = getr2(d,verbose);

% and parse back fields
if size(timecourses,1)==1
  outd.r2 = squeeze(outd.r2(1,1,1,:,:));
  outd.ehdr = squeeze(outd.ehdr(1,1,1,:,:));
  outd.ehdrste = squeeze(outd.ehdrste(1,1,1,:,:));
else
  outd.r2 = squeeze(outd.r2(1,1,:,:,:));
  outd.ehdr = squeeze(outd.ehdr(1,1,:,:,:));
  outd.ehdrste = squeeze(outd.ehdrste(1,1,:,:,:));
end

% get dimensions right
if nhdr == 1
  outd.ehdr = outd.ehdr';
  outd.ehdrste = outd.ehdrste';
end

if ~ieNotDefined('tr')
  outd.time = tr/2:tr:(hdrlen*tr);
end

outd.covar = pinv(d.scm'*d.scm);
outd = rmfield(outd,'data');


