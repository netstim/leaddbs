% eventRelatedTSeries.m
%
%        $Id$
%      usage: eventRelatedTSeries.m(timecourses,nhdr,hdrlen,scm)
%         by: justin gardner
%       date: 05/25/07
%    purpose: function that will use getr2 (used in mrLoadRet) to do the event related
%             analysis on a time series or an array of time series
%             e.g. to get the event related analysis for a mena ROI
%             time series
%
%             v = newView;
%             d = erAnal2d(1,3);
%             roi = loadROITSeries(v,'l_v1',1,3);
%             tSeries = mean(roi.tSeries);
%             d = eventRelatedTSeries(tSeries,d.nhdr,d.hdrlen,d.scm)
%             deleteView(v);
%
function outd = eventRelatedTSeries(timecourses,nhdr,hdrlen,scm)

% check arguments
if ~any(nargin == [4])
  help eventRelatedTSeries
  return
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
outd = getr2(d);

% and parse back fields
if size(timecourses,1)
  outd.r2 = squeeze(outd.r2(1,1,1,:));
  outd.ehdr = squeeze(outd.ehdr(1,1,1,:,:));
  outd.ehdrste = squeeze(outd.ehdrste(1,1,1,:,:));
else
  outd.r2 = squeeze(outd.r2(1,1,:,:,:));
  outd.ehdr = squeeze(outd.ehdr(1,1,:,:,:));
  outd.ehdrste = squeeze(outd.ehdrste(1,1,:,:,:));
end
